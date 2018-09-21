"""
	struct DCAgraph
"""
struct DCAgraph
    J::Array{Float64,2}
    h::Array{Float64,1}
    L::Int64
    q::Int64
end

"""
	function readparam(infile::String ; format="mcmc", q=0)

Read dca parameters from `infile`. Format option can be either 
- `"mcmc"`: `J i j a b val`
- `"mat"`: One line of `infile` represents the vector `J[i,a][:]`. This is useful for parameters stored in dlm format. Optional argument `q` is needed in this case. 
Output is of type `DCAgraph`.
"""
function readparam(infile::String ; format="mcmc", q=0)
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
	function readparammcmc(infile::String)

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
	hlines = lines[find(m->m!=nothing,map(x->match(r"^h",x),lines))] # step 1
	L = findmax(map(l->parse(split(l," ")[2]),hlines))[1]+1 # step 2 for `L`
	q = findmax(map(l->parse(split(l," ")[3]),hlines))[1]+1 # step 2 for `q`

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
        x = convert(Array{Int64,1},x[1:end-1]) + 1
            if in('h',l)
                h[(x[1]-1)*q + x[2]] = v
            elseif in('J',l)
                J[(x[1]-1)*q + x[3], (x[2]-1)*q + x[4]] = v
                J[(x[2]-1)*q + x[4], (x[1]-1)*q + x[3]] = v
            end
	end 
	return (J,h,L,q)
end