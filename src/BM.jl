module BM

using DCATools
using DCATools.MutLand
using Printf
using Statistics
using Parameters

export BMmeta, bmlearn

include("objects.jl")
include("BM/computegradient.jl")
include("BM/makebm.jl")

end
