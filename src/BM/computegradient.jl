"""
	updateparameters!(g::DCAGraph, grad::DCAGrad)

Add gradient `grad` to graph `g`, modifying it. 
"""
function updateparameters!(g::DCAGraph, grad::DCAGrad)
	g.J .+= grad.stepJ .* grad.gradJ
	g.h .+= grad.steph .* grad.gradh
end


"""
	computegradient!(grad::DCAGrad, sample::Array{Int64,2}, f1::Array{Float64,1}, f2::Array{Float64,2})

Compute gradient corresponding to difference between pairwise frequencies measured in `sample` and targets `f1` and `f2`. Modify input `grad` in order to avoid allocation. Step size are set to 0.
"""
function computegradient!(grad::DCAGrad, sample::Array{Int64,2}, f1::Array{Float64,1}, f2::Array{Float64,2})
	p1, p2 = pairwise_frequencies(sample)
	grad.gradJ = f2-p2
	grad.gradh = f1-p1
	grad.stepJ .*= 0
	grad.steph .*= 0
	return grad,p1,p2
end

"""
	computegradient(sample::Array{Int64,2}, f1::Array{Float64,1}, f2::Array{Float64,2}, q::Int64)

Compute gradient corresponding to difference between pairwise frequencies measured in `sample` and targets `f1` and `f2`. Allocate a new `DCAGrad` object. Step size are set to 0.
"""
function computegradient(sample::Array{Int64,2}, f1::Array{Float64,1}, f2::Array{Float64,2}, q::Int64)
	p1, p2 = pairwise_frequencies(sample)
	L = size(sample,2)

	grad = DCAGrad(L,q)
	grad.gradJ = f2-p2
	grad.gradh = f1-p1

	return grad, p1, p2
end

"""
	computegradient(md::MutData, mapping::Dict{Float64, Float64}, meta::BMmeta)

Gradient due to differences between measured fitness and energies. Compute differences between energies and fitness over all mutants in `md`. Result is scaled by sample size and by λ. 
"""
function computegradient(md::MutData, mapping::Dict{Float64, Float64}, meta::BMmeta)
	grad = DCAGrad(md.L,md.q)

	for mut in md.mutant
		if size(mut.smut,1) == 1 # gradient on fields
			i = mut.smut[1].i
			a = mut.smut[1].a
			grad.gradh[(i-1)*md.q + a] += mut.fitness - mapping[mut.E]
		else size(mut.smut,1) == 2 # gradient on coupling
			i = mut.smut[1].i
			a = mut.smut[1].a
			j = mut.smut[2].j
			b = mut.smut[2].b
			grad.gradJ[(i-1)*md.q + a, (j-1)*md.q + b] += mut.fitness - mapping[mut.E] 
		end
	end
	grad.gradJ = grad.gradJ + grad.gradJ'
	grad.gradh .*= meta.integrative_lambda / (1 - meta.integrative_lambda) / meta.integrative_M
	grad.gradJ .*= meta.integrative_lambda / (1 - meta.integrative_lambda) / meta.integrative_M
	return grad
end


"""
	computel2(g::DCAGraph, lambda::Float64)

Return gradient corresponding to l2 regularization. 
"""
function computel2(g::DCAGraph, lambda::Float64)

	grad = DCAGrad(g.L, g.q)
	grad.gradJ = -lambda*g.J
	grad.gradh = -lambda*g.h

	return grad
end
"""
	computel2(g::DCAGraph, lambdaJ::Float64, lambdah::Float64)

Return gradient corresponding to l2 regularization. 
"""
function computel2(g::DCAGraph, lambdaJ::Float64, lambdah::Float64)

	grad = DCAGrad(g.L, g.q)
	grad.gradJ = -lambdaJ*g.J
	grad.gradh = -lambdah*g.h

	return grad
end

"""
	computel1(g::DCAGraph, lambda::Float64)

Return gradient corresponding to l1 regularization. Experimental, should not be used. 
"""
function computel1!(g::DCAGraph, cgrad::DCAGrad, lambda::Float64)
	# If J = 0 and the gradient is smaller than lambda, then J should not move 
	# Else, if J != 0, add lambda to the gradient
	grad = DCAGrad(g.L, g.q)
	eps = lambda
	for j in 1:g.L
		for i in 1:g.L
			if i!=j
				for b in 1:g.q
					for a in 1:g.q
						if g.J[(i-1)*g.q + a, (j-1)*g.q + b]>eps
							grad.gradJ[(i-1)*g.q + a, (j-1)*g.q + b] =  -min(lambda, g.J[(i-1)*g.q + a, (j-1)*g.q + b])
						elseif g.J[(i-1)*g.q + a, (j-1)*g.q + b]<eps
							grad.gradJ[(i-1)*g.q + a, (j-1)*g.q + b] =  min(lambda, -g.J[(i-1)*g.q + a, (j-1)*g.q + b])
						else
							if abs(cgrad[(i-1)*g.q + a, (j-1)*g.q + b]) < eps
								cgrad[(i-1)*g.q + a, (j-1)*g.q + b] = 0
							end
						end
					end
				end
			end
		end
	end

	# grad.gradJ = -lambda*(g.J .> 0) + lambda*(g.J.<0)
	# grad.gradh = -lambda*(g.h .> 0) + lambda*(g.h.<0)

	return grad
end


"""
	computestepsize!(newgrad::DCAGrad, prevgrad::DCAGrad, meta::BMmeta)

Update step size of the gradient descent. For single parameters (i,j,a,b), the difference of sign between the current gradient and the previous one is used. If both are of the same sign, step size is increased. Other wise, it is decreased. 
"""
function computestepsize!(newgrad::DCAGrad, prevgrad::DCAGrad, meta::BMmeta)
	xJ = (newgrad.gradJ .* prevgrad.gradJ) .>=0
	newgrad.stepJ .= prevgrad.stepJ .* ((.!xJ * meta.adaptsJdown) + (xJ * meta.adaptsJup))
	newgrad.stepJ .= min.(newgrad.stepJ, meta.stepJmax)

	xh = (newgrad.gradh .* prevgrad.gradh) .>=0
	newgrad.steph .= prevgrad.steph .* ((.!xh * meta.adaptshdown) + (xh * meta.adaptshup))
	newgrad.steph .= min.(newgrad.steph, meta.stephmax)
end

"""
"""
function gradnorm(grad::DCAGrad)
	return sum(grad.gradJ.^2) + sum(grad.gradh.^2), sum(grad.gradh.^2), sum(grad.gradJ.^2)
end

"""
	gradconst(prevgrad::DCAGrad, newgrad::DCAGrad)

Scalar product between two gradients, normalized by their norms. Defines cosine between the two vectors. 
"""
function gradconst(prevgrad::DCAGrad, newgrad::DCAGrad)
	outJ = sum(prevgrad.gradJ .* newgrad.gradJ) / sqrt(sum(newgrad.gradJ.^2) * sum(prevgrad.gradJ.^2))
	outh = sum(prevgrad.gradh .* newgrad.gradh) / sqrt(sum(newgrad.gradh.^2) * sum(prevgrad.gradh.^2))
	return outh,outJ
end
