export computegradient, computel2


"""
	updateparameters!(g::DCAgraph, grad::DCAgrad)

Add gradient `grad` to graph `g`, modifying it. 
"""
function updateparameters!(g::DCAgraph, grad::DCAgrad)
	g.J .+= grad.stepJ .* grad.gradJ
	g.h .+= grad.steph .* grad.gradh
end


"""
	computegradient!(grad::DCAgrad, sample::Array{Int64,2}, f1::Array{Float64,1}, f2::Array{Float64,2})

Compute gradient corresponding to difference between pairwise frequencies measured in `sample` and targets `f1` and `f2`. Modify input `grad` in order to avoid allocation. Step size are set to 0.
"""
function computegradient!(grad::DCAgrad, sample::Array{Int64,2}, f1::Array{Float64,1}, f2::Array{Float64,2})
	p1, p2 = computefreqs(sample)
	grad.gradJ = f2-p2
	grad.gradh = f1-p1
	grad.stepJ .*= 0
	grad.steph .*= 0
	return grad,p1,p2
end

"""
	computegradient(sample::Array{Int64,2}, f1::Array{Float64,1}, f2::Array{Float64,2}, q::Int64)

Compute gradient corresponding to difference between pairwise frequencies measured in `sample` and targets `f1` and `f2`. Allocate a new `DCAgrad` object. Step size are set to 0.
"""
function computegradient(sample::Array{Int64,2}, f1::Array{Float64,1}, f2::Array{Float64,2}, q::Int64)
	p1, p2 = computefreqs(sample)
	L = size(sample,2)

	grad = DCAgrad(L,q)
	grad.gradJ = f2-p2
	grad.gradh = f1-p1

	return grad, p1, p2
end


"""
	computel2!(grad::DCAgrad, g::DCAgraph, lambda::Float64)
"""
function computel2!(grad::DCAgrad, g::DCAgraph, lambda::Float64)
	grad.gradJ = -lambda*g.J
	grad.gradh = -lambda*g.h
	grad.stepJ .*= 0
	grad.steph .*= 0
	return grad
end

"""
	computel2(g::DCAgraph, lambda::Float64)

Return gradient corresponding to l2 regularization. 
"""
function computel2(g::DCAgraph, lambda::Float64)

	grad = DCAgrad(g.L, g.q)
	grad.gradJ = -lambda*g.J
	grad.gradh = -lambda*g.h

	return grad
end

"""
	computel1(g::DCAgraph, lambda::Float64)

Return gradient corresponding to l1 regularization. 
"""
function computel1(g::DCAgraph, lambda::Float64)

	grad = DCAgrad(g.L, g.q)
	grad.gradJ = -lambda*(g.J .> 0) + lambda*(g.J.<0)
	# grad.gradh = -lambda*(g.h .> 0) + lambda*(g.h.<0)
	# display(grad.gradJ)

	return grad
end


"""
	computestepsize!(newgrad::DCAgrad, prevgrad::DCAgrad, meta::BMmeta)

Update step size of the gradient descent. For single parameters (i,j,a,b), the difference of sign between the current gradient and the previous one is used. If both are of the same sign, step size is increased. Other wise, it is decreased. 
"""
function computestepsize!(newgrad::DCAgrad, prevgrad::DCAgrad, meta::BMmeta)
	xJ = (newgrad.gradJ .* prevgrad.gradJ) .>=0
	newgrad.stepJ .= prevgrad.stepJ .* ((.!xJ * meta.adaptsJdown) + (xJ * meta.adaptsJup))
	newgrad.stepJ .= min.(newgrad.stepJ, meta.stepJmax)

	xh = (newgrad.gradh .* prevgrad.gradh) .>=0
	newgrad.steph .= prevgrad.steph .* ((.!xh * meta.adaptshdown) + (xh * meta.adaptshup))
	newgrad.steph .= min.(newgrad.steph, meta.stephmax)
end

"""
"""
function gradnorm(grad::DCAgrad)
	return sum(grad.gradJ.^2) + sum(grad.gradh.^2), sum(grad.gradh.^2), sum(grad.gradJ.^2)
end

"""
	gradconst(prevgrad::DCAgrad, newgrad::DCAgrad)

Scalar product between two gradients, normalized by their norms. Defines cosine between the two vectors. 
"""
function gradconst(prevgrad::DCAgrad, newgrad::DCAgrad)
	outJ = sum(prevgrad.gradJ .* newgrad.gradJ) / sqrt(sum(newgrad.gradJ.^2) * sum(prevgrad.gradJ.^2))
	outh = sum(prevgrad.gradh .* newgrad.gradh) / sqrt(sum(newgrad.gradh.^2) * sum(prevgrad.gradh.^2))
	return outh,outJ
end
