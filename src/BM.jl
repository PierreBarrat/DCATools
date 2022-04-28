module BM

using DCATools
using DCATools.MutLand
using DCATools.MCMC
using Printf
using Statistics
using Parameters

## Things I could improve
#=
- The update_tau thing is a bit weird
- Check that old scripts work with this version
=#

export DCAgrad, gradequal, BMmeta

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
    *(λ::Float64, g::DCAgrad)
"""
function *(λ::Real, g::DCAgrad)
    return DCAgrad(λ * g.gradJ, λ*g.gradh, g.stepJ, g.steph, g.L, g.q)
end
"""
    *(g::DCAgrad, λ::Real)
"""
function *(g::DCAgrad, λ::Float64)
    return DCAgrad(λ * g.gradJ, λ*g.gradh, g.stepJ, g.steph, g.L, g.q)
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

Meta parameters for the BM learning. Those are never modified in the course of the computation. Call `BMmeta(;kwargs...)` for constructing an instance of `BMmeta`.

### Regularization
- `l2J::Float64 = 0.01`: l2 regularization for couplings
- `l2h::Float64 = 0.01`: l2 regularization for fields
- `l2::Union{Missing, Float64} = missing`: joint l2 regularization. Supersedes `l2J` and `l2h`. 
- `l1::Float64 = 0.`: l1 regularization
### Estimating gradient
- `Minit::Int64 = 1_000`: initial size of MCMC sample
- `Mmax::Int64 = 100_000`: maximal size of MCMC samples to estimate gradient
- `adaptMup::Float64 = 1.05`: factor to increase sample sizes when gradient fluctuates too much. 
### Gradient ascent step
- `basestepJ::Float64 = 0.05`: initial stepsize for couplings
- `basesteph::Float64 = 0.05`: *idem* for fields
- `stepJmax::Float64 = 0.5`: max. stepsize for couplings
- `stephmax::Float64 = 0.5`: *idem* for fields
- `adaptsJup::Float64 = 1.2`: factor to increase stepsize for couplings
- `adaptsJdown::Float64 = 0.5`: factor to decrease stepsize for couplings
- `adaptshup::Float64 = 1.2`: *idem* for fields
- `adaptshdown::Float64 = 0.5`: *idem* for fields
### Integrating mutational data
- `integrative_lambda::Float64 = 0.`: Strength of integration, between `0` and `1`. `0` ignores mutational data. 
- `integrative_M::Int64 = 1`: Size of alignment that target frequencies were computed with. Not used if no mutational data. 
### Saving results
- `saveparam::Int64 = 10`: number of iterations between parameter save. 
- `logfile::String = "log.txt"`: store log of bmlearn
- `infofile::String = "info.txt"`: store information about bmlearn run
- `comment::String = ""`: comment to be stored in `infofile`. 
- `savefolder::String = "./"`: created at runtime. Set to `""` for no file output.
### Misc. 
- `nit::Int64 = 50`: number of iterations
- `nprocs::Int64 = 1`: parallel MCMC is used if `nprocs>1` (not working super well if I remember correctly).
- `update_tau::Int64 = 10`: steps separating MCMC samples is re-estimated every `update_tau`. 
- `verbose::Bool = false`: verbosity
"""
@with_kw struct BMmeta
    # Regularization
    l2J::Float64 = 0.01
    l2h::Float64 = 0.01
    l2::Union{Missing, Float64} = missing
    l1::Float64 = 0.
    # Estimating gradient
    Minit::Int64 = 1_000
    Mmax::Int64 = 100_000
    adaptMup::Float64 = 1.05
    # Gradient ascent step
    basestepJ::Float64 = 0.05
    basesteph::Float64 = 0.05
    stepJmax::Float64 = 0.5
    stephmax::Float64 = 0.5
    adaptsJup::Float64 = 1.2
    adaptsJdown::Float64 = 0.5
    adaptshup::Float64 = 1.2
    adaptshdown::Float64 = 0.5
    # Integrating mutational data
    integrative_lambda::Float64 = 0.
    integrative_M::Int64 = 1
    # Saving results
    saveparam::Int64 = 10
    logfile::String = "log.txt"
    infofile::String = "info.txt"
    comment::String = ""
    savefolder::String = "./"
    # Misc. 
    nit::Int64 = 50
    nprocs::Int64 = 1
    update_tau::Int64 = 10
    verbose::Bool = false
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


include("BM/computegradient.jl")
include("BM/makebm.jl")

end
