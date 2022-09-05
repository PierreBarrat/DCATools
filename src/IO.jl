"""
	DCAGraph(infile::AbstractString, format=:extended; q=0)

Read dca parameters from `infile` into a `DCAGraph` object.
Format option can be either:
- `:extended`: `J i j a b val`
- `:matrix`: One line of `infile` represents the vector `J[i,a][:]`. This is useful for parameters stored in dlm format. Optional argument `q` is needed in this case.
"""
function DCAGraph(infile::AbstractString, format=:extended; q=0)
	g = if format == :matrix
		@assert q > 0 "With `:mat` option, you need to specify `q`. Received q=$q"
		J, h, L = read_graph_matrix(infile, q)
		DCAGraph(J, h, L, q)
	elseif format == :extended
		read_graph_extended(infile)
	else
		error("Unrecognized format argument. Options are `:extended` or `:mat`.\n")
	end
	return g
end


"""
    read_graph_matrix(file::AbstractString, q::Int)

Reads potts parameters in dlm format (*ie* stored as a matrix, with h being the last line).  
Output `J`,`h` and `L`. 
"""
function read_graph_matrix(file::AbstractString, q::Int)
    f = open(file)
    try 
    	J::Array{Float64,2} = readdlm(f,Float64)
    catch err
    	println("error when attempting to dlmread file $file\n")
    	error("$err")
    end
    close(f)
    h::Array{Float64,1} = J[end,:]
    J = J[1:(end-1),:]

    # Checking for size 
    if size(h,1)%q != 0
    	error("Incorrect file size ; q=$q, number of cols=$(size(h,1))\n")
	else
		L = Int(size(h,1)/q)
	end
    return (J,h,L)
end

function read_graph_extended(file)
	## Go through file twice: first to get L and q, the second to store parameters
	format = 1
	q = 0
	L = 0
	for line in eachline(file)
		@assert is_valid_line(line) "Format problem with line:\n $line \n--> Expected `J i j a b` or `h i a`."
		if !isempty(line) && line[1] == 'h'
			i, a, val = parse_field_line(line)
			if i > L
				L = i
			end
			if a > q
				q = a
			end
			if i == 0 || q == 0
				format = 0
			end
		end
	end
	if format == 0
		L += 1
		q += 1
	end

	J = zeros(Float64, L*q, L*q)
	h = zeros(Float64, L*q)
	g = DCAGraph(J, h, L, q)
	for line in eachline(file)
		if line[1] == 'J'
			i, j, a, b, val = parse_coupling_line(line)
			format == 1 ? (g[i, j, a, b] = val) : (g[i+1, j+1, a+1, b+1] = val)
		elseif line[1] == 'h'
			i, a, val = parse_field_line(line)
			format == 1 ? (g[i, a] = val) : (g[i+1, a+1] = val)
		end
	end

	return g
end

function is_valid_line(line)
	if isnothing(match(r"J [0-9]+ [0-9]+ [0-9]+ [0-9]+", line)) &&
		isnothing(match(r"h [0-9]+ [0-9]+", line))
		return false
	else
		return true
	end
end
function parse_field_line(line)
	s = split(line, " ")
	i = parse(Int, s[2])
	a = parse(Int, s[3])
	val = parse(Float64, s[4])
	return i, a, val
end
function parse_coupling_line(line)
	s = split(line, " ")
	i = parse(Int, s[2])
	j = parse(Int, s[3])
	a = parse(Int, s[4])
	b = parse(Int, s[5])
	val = parse(Float64, s[6])
	return i, j, a, b, val
end


"""
	writeparam(outfile::AbstractString, g::DCAGraph; format="mat")

Write graph `g` to file `outfile`: 
- as a matrix if `format=="mat"`
- using `J i j a b value` if `format=="mcmc"`
"""
function writeparam(outfile::AbstractString, g::DCAGraph; format="mat")
	if format == "mat"
		writedlm(outfile, [round.(g.J, digits = 5) ; round.(g.h', digits = 5)], " ")
	elseif format == "mcmc"
		writeparammcmc(outfile,g)
	else
		error("inputoutput.jl - writeparam: Unrecognized format argument.\n")
	end
end


"""
	writeparammcmc(outfile::AbstractString, g::DCAGraph)

Write graph `g` to file `outfile` using format `J i j a b value`.
"""
function writeparammcmc(outfile::AbstractString, g::DCAGraph)
	f = open(outfile, "w")
	for i in 1:g.L
		for j in (i+1):g.L
			for a in 1:g.q
				for b in 1:g.q
					write(f, "J $(i-1) $(j-1) $(a-1) $(b-1) $(g.J[(j-1)*g.q+b, (i-1)*g.q+a])\n" )
				end
			end
		end
	end
	for i in 1:g.L
		for a in 1:g.q
			write(f, "h $(i-1) $(a-1) $(g.h[(i-1)*g.q+a])\n" )
		end
	end
	close(f)
end



"""
	readmsanum(infile::AbstractString ; format=1, header=false)

Read an MSA stored in `infile` in a numerical format. 

If `format=1`, amino acids should be mapped from 1 to `q`. If `format=0`, they should be mapped from 0 to `q-1`.
`header` argument allows for discarding the first line of `infile`. 
"""
function readmsanum(infile::AbstractString ; format=1, header=false)
	Y = Array{Float64,2}(undef,0,0)
	try 
		if header
			Y = readdlm(infile, Int, skipstart=1)
		else
			Y = readdlm(infile, Int)
		end	
	catch err
		println("inputoutput.jl - readmsanum: readdlm failed on file $infile. The alignment may not be of the correct format.")
		error(err)
	end

	if format==0
		Y .+= 1
	elseif format==1
		if findmin(Y)[1] == 0
			error("inputoutput.jl - readmsanum: file $infile contains a 0 value. `format` should be set to 0.")
		end
	end

	return Y
end
