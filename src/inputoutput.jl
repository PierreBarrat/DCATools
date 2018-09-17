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
	if format=="mat"
		if q==0
			error("inputoutput.jl - readparam: With \"mat\" option, you need to specify `q`.\n")
		else
			(J,h,L) = readparammat(infile,q)
		end
	else
		error("inputoutput.jl - readparam: Unrecognized format argument.\n")
	end

	return DCAgraph(J,h,L,q)
end


"""
    readparammat(infile::String, q::Int64)

Reads potts parameters in dlm format (ie stored as a matrix, with h being the last line).  
Outputs J and h. 
"""
function readparammat(infile::String, q::Int64)
    f = open(parameters)
    J::Array{Float64,2} = readdlm(f,Float64)
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