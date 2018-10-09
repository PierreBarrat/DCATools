###


export readparam, writeparam, readmsanum

###
"""
	readparam(infile::String ; format="mat", q=0)

Read dca parameters from `infile`. Format option can be either 
- `"mcmc"`: `J i j a b val`
- `"mat"`: One line of `infile` represents the vector `J[i,a][:]`. This is useful for parameters stored in dlm format. Optional argument `q` is needed in this case. 
Output is of type `DCAgraph`.
"""
function readparam(infile::String ; format="mat", q=0)
	if format == "mat"
		if q==0
			error("inputoutput.jl - readparam: With \"mat\" option, you need to specify `q`.\n")
		else
			(J,h,L) = readparammat(infile,q)
		end
	elseif format == "mcmc"
		(J,h,L,q) = readparammcmc(infile)
	else
		error("inputoutput.jl - readparam: Unrecognized format argument.\n")
	end
	return DCAgraph(J,h,L,q)
end


"""
    readparammat(infile::String, q::Int64)

Reads potts parameters in dlm format (*ie* stored as a matrix, with h being the last line).  
Output `J`,`h` and `L`. 
"""
function readparammat(infile::String, q::Int64)
    f = open(infile)
    try 
    	J::Array{Float64,2} = readdlm(f,Float64)
    catch error
    	error("inputoutput.jl - readparammat: error when attempting to dlmread file $infile\n")
    end
    close(f)
    h::Array{Float64,1} = J[end,:]
    J = J[1:(end-1),:]

    # Checking for size 
    if size(h,1)%q != 0
    	error("inputoutput.jl - readparammat: Incorrect file size ; q=$q, number of cols=$(size(h,1))\n")
	else
		L = Int64(size(h,1)/q)
	end
    return (J,h,L)
end

"""
	readparammcmc(infile::String)

Read potts parameters in mcmc format: `J i j a b value` or `h i a value`. 
Output `J`, `h` and `L`.
"""
function readparammcmc(infile::String)
	lines = open(infile) do f 
		lines = readlines(f)
	end	

	## Strategy
	# 1. Get lines starting with h
	# 2. Retrieve value for `L` and `q` with those. 
	hlines = lines[findall(m->m!=nothing,map(x->match(r"^h",x),lines))] # step 1
	L = findmax(map(l->Meta.parse(split(l," ")[2]),hlines))[1]+1 # step 2 for `L`
	q = findmax(map(l->Meta.parse(split(l," ")[3]),hlines))[1]+1 # step 2 for `q`

	## Now we can attribute memory
	J = zeros(Float64, L*q, L*q)
	h = zeros(Float64, L*q)
	# i = 0
	# j = 0
	# a = 0
	# b = 0
	for l in lines
		x = readdlm(IOBuffer(l[2:end]))
        v = x[end]
        x = convert(Array{Int64,1},x[1:end-1]) .+ 1
            if in('h',l)
                h[(x[1]-1)*q + x[2]] = v
            elseif in('J',l)
                J[(x[1]-1)*q + x[3], (x[2]-1)*q + x[4]] = v
                J[(x[2]-1)*q + x[4], (x[1]-1)*q + x[3]] = v
            end
	end 
	return (J,h,L,q)
end



"""
	writeparam(outfile::String, g::DCAgraph; format="mat")

Write graph `g` to file `outfile`: 
- as a matrix if `format=="mat"`
- using `J i j a b value` if `format=="mcmc"`
"""
function writeparam(outfile::String, g::DCAgraph; format="mat")
	if format == "mat"
		writedlm(outfile, [g.J ; g.h'], " ")
	elseif format == "mcmc"
		writeparammcmc(outfile,g)
	else
		error("inputoutput.jl - writeparam: Unrecognized format argument.\n")
	end
end


"""
	writeparammcmc(outfile::String, g::DCAgraph)

Write graph `g` to file `outfile` using format `J i j a b value`.
"""
function writeparammcmc(outfile::String, g::DCAgraph)
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
	readmsanum(infile::String ; format=1, header=false)

Read an MSA stored in `infile` in a numerical format. 

If `format=1`, amino acids should be mapped from 1 to `q`. If `format=0`, they should be mapped from 0 to `q-1`.
`header` argument allows for discarding the first line of `infile`. 
"""
function readmsanum(infile::String ; format=1, header=false)
	Y = Array{Float64,2}(undef,0,0)
	try 
		if header
			Y = readdlm(infile, Int64, skipstart=1)
		else
			Y = readdlm(infile, Int64)
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