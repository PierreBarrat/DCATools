module DCABM

using DCATools
using DCAMutLand
using DCAMCMC
using Printf
using Statistics

export DCAgrad, gradequal

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
function DCAgrad(L::Int64, q::Int64)
    return DCAgrad(zeros(Float64, L*q,L*q), zeros(Float64,L*q), zeros(Float64, L*q,L*q), zeros(Float64,L*q), L, q)
end

import Base: +, *


"""
	+(A::DCAgrad, B::DCAgrad)

Add two gradients. Step size is the one of the first argument. 
"""
function +(A::DCAgrad, B::DCAgrad)
    return DCAgrad(A.gradJ + B.gradJ, A.gradh + B.gradh, A.stepJ, A.steph, A.L, A.q)
end

"""
    gradequal(x, y)

Test equality between two `DCAgrad` objects.
"""
function gradequal(x, y)
    return (x.gradJ==y.gradJ) && (x.gradh==y.gradh) && (x.stepJ==y.stepJ) && (x.steph==y.steph)
end


"""
    BMmeta

Meta parameters for the BM learning. Those are never modified in the course of the computation. 
"""
struct BMmeta
    # Regularization
    l2::Float64
    l1::Float64
    # Step size
    basestepJ::Float64
    basesteph::Float64
    stepJmax::Float64
    stephmax::Float64
    adaptsJup::Float64
    adaptsJdown::Float64
    adaptshup::Float64
    adaptshdown::Float64
    # 
    adaptMup::Float64
    Mmax::Int64
    #
    lambda::Float64
    Mmsa::Int64
    #
    saveparam::Int64
end

"""
    BMlog

Log of BM learning. Different informations about the current state of the learning process are stored here. 
"""
mutable struct BMlog
    samplesize
    gradnorm
    gradnormh
    gradnormJ
    gradconsth
    gradconstJ
    tau
    corcor
    slopecor
    cormag
    cormutants
end
function BMlog()
    return BMlog(0,0,0,0,0,0,0,0,0,0,0)
end


include("DCABM/computegradient.jl")
include("DCABM/makebm.jl")

end