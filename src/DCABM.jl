module DCABM

using DCATools
using DCAMutland

export DCAgrad, updateparameters!

"""
	DCAgrad

Store direction of gradient in `dJ` and `dh`, and step size in `stepJ` and `steph`. 
"""
mutable struct DCAgrad
    gradJ::Array{Float64,2}
    gradh::Array{Float64,1}
    stepJ::Array{Float64,2}
    steph::Array{Float64,1}
    L::Int64
    q::Int64
end

import Base: +

"""
	+(A::DCAgrad, B::DCAgrad)

Add two gradients. Step size is the one of the first argument. 
"""
function +(A::DCAgrad, B::DCAgrad)
    return DCAgrad(A.gradJ + B.gradJ, A.gradh + B.gradh, A.stepJ, A.steph, A.L, A.q)
end

"""
	updateparameters!(g::DCAgraph, grad::DCAgrad)

Add gradient `grad` to graph `g`, modifying it. 
"""
function updateparameters!(g::DCAgraph, grad::DCAgrad)
	g.J .+= grad.stepJ .* grad.gradJ
	g.h .+= grad.steph .* grad.gradh
end

include("DCABM/computegradient.jl")

end