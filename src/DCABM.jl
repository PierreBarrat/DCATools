module DCABM

using DCATools
using DCAMutland

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

"""
	updateparameters!(g::DCAgraph, grad::DCAgrad)

Add gradient `grad` to graph `g`, modifying it. 
"""
function updateparameters!(g::DCAgraph, grad::DCAgrad)
	g.J .+= grad.stepJ .* grad.gradJ
	g.h .+= grad.steph .* grad.gradh
end

end