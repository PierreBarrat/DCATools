"""
	computegradient(sample::Array{Int64,2}, f1::Array{Float61,1}, f2::Array{Float64,2})

Compute gradient corresponding to difference between pairwise frequencies measured in `sample` and targets `f1` and `f2`. Modify input `grad` in order to avoid allocation. Step size are set to 0.
"""
function computegradient!(grad::DCAgrad, sample::Array{Int64,2}, f1::Array{Float61,1}, f2::Array{Float64,2})
	p1, p2 = computefreqs(sample)
	grad.gradJ = f2-p2
	grad.gradh = f1-p1
	grad.stepJ .*= 0
	grad.steph .*= 0
	return grad
end

"""
	computel2(g::DCAgraph, lambda::Float64)
"""
function computel2!(grad::DCAgrad, g::DCAgraph, lambda::Float64)
	grad.gradJ = -lambda*g.J
	grad.gradh = -lambda*g.h
	grad.stepJ .*= 0
	grad.steph .*= 0
	return grad
end

