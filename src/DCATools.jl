module DCATools

using DelimitedFiles
using LinearAlgebra
using Printf
using Random
using Statistics

import Base: *, getindex, setindex!, size


include("objects.jl")
export DCAgraph

include("global.jl")
include("inputoutput.jl")
export writeparam, readmsanum

include("alignment_frequencies.jl")
export pairwise_frequencies, computeweights

include("modeltools.jl")
export switchgauge!, computeenergies

include("contactprediction.jl")
export PPV, Fapc

include("misc.jl")
export fitquality, threepointscor, corr3p

include("MF.jl")

# Other sub-modules
include("sampling.jl")
export sample


#
include("MutLand/MutLand.jl")
export MutData
export Mutant
export SingleMut

export computeenergies!
export computeepistasis!
export findsinglemut
export finddoublemut
export mapenergies!
export mapsinglemut!
export readmutdata

#
include("BM.jl")
export bmlearn


end
