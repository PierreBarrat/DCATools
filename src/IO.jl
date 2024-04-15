"""
	DCAGraph(infile::AbstractString, format=:extended; q=0)

Read dca parameters from `infile` into a `DCAGraph` object.
Format option can be either:
- `:extended`: `J i j a b val`
- `:matrix`: One line of `infile` represents the vector `J[i,a][:]`.
	This is useful for parameters stored in dlm format. Optional argument `q` is needed in this case.
"""
function DCAGraph(infile::AbstractString, format=:extended; q=0)
	g = if format == :matrix
		@assert q > 0 "With `:matrix` option, you need to specify `q`. Got q=$q"
		read_graph_matrix(infile, q)
	elseif format == :extended
		read_graph_extended(infile)
	else
		error("Unrecognized format argument. Options are `:extended` or `:matrix`.\n")
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
    return DCAGraph(;J, h, L, q)
end

function read_graph_extended(file)
	## Go through file twice: first to get L and q, the second to store parameters
	q = 0
	L = 0
	min_idx = Inf
	index_style = 1
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
			if i < min_idx || a < min_idx
				min_idx = min(i, a)
			end
			# if i == 0 || q == 0
			# 	index_style = 0
			# end
		end
	end
	@assert min_idx == 1 || min_idx == 0 "Issue with indexing: smallest index found is $min_idx"
	index_style = (min_idx == 0 ? 0 : 1)
	if index_style == 0
		L += 1
		q += 1
	end

	g = DCAGraph(;L, q)
	for line in eachline(file)
		if line[1] == 'J'
			i, j, a, b, val = parse_coupling_line(line)
			index_style == 0 && (i += 1; j += 1; a += 1; b += 1)
			g[i, j, a, b] = val
			g[j, i, b, a] = val
		elseif line[1] == 'h'
			i, a, val = parse_field_line(line)
			index_style == 0 && (i += 1; a += 1)
			g[i, a] = val
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
	write(file::AbstractString, g::DCAGraph, format = :extended; sigdigits=5)

Write parameters of `g` to `file`:
- using `J i j a b value` if `format == :extended`
- as a matrix if `format == :matrix`
"""
function write(
	file::AbstractString, g::DCAGraph, format = :extended;
	sigdigits=5, index_style=1,
)
	if format == :matrix
		writedlm(
			file,
			[round.(g.J, digits = sigdigits) ; round.(g.h', digits = sigdigits)],
			" ",
		)
	elseif format == :extended
		write_graph_extended(file, g, sigdigits, index_style)
	else
		error("inputoutput.jl - writeparam: Unrecognized format argument.\n")
	end
end


"""
	write_graph_extended(file::AbstractString, g::DCAGraph; sigdigits=5)

Write graph `g` to file `file` using format `J i j a b value`.
"""
function write_graph_extended(file::AbstractString, g::DCAGraph, sigdigits, index_style)
	@assert index_style == 0 || index_style == 1 "Got `index_style==`$(index_style)"
	f = open(file, "w")
	for i in 1:g.L
		for j in (i+1):g.L
			for a in 1:g.q
				for b in 1:g.q
					val = round(g[j,i,b,a]; sigdigits)
					if index_style == 0
						write(f, "J $(i-1) $(j-1) $(a-1) $(b-1) $val\n")
					elseif index_style == 1
						write(f, "J $i $j $a $b $val\n")
					end
				end
			end
		end
	end
	for i in 1:g.L
		for a in 1:g.q
			val = round(g[i,a]; sigdigits)
			if index_style == 0
				write(f, "h $(i-1) $(a-1) $val\n" )
			elseif index_style == 1
				write(f, "h $i $a $val\n" )
			end
		end
	end
	close(f)
end

function write(file::AbstractString, model::ProfileModel; sigdigits=5, delim='\t')
    L, q = model.L, model.q
    header = reshape([model.mapping[a] for a in 1:q], 1, q)
    P = round.(reshape(model.w, q, L); sigdigits)'
    writedlm(file, [header; P], delim)
end




