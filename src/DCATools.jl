module DCATools

using DelimitedFiles
using FASTX
using LinearAlgebra
using Printf
using Random
using Statistics

import Base: *, firstindex, lastindex, getindex, iterate, length, eltype, setindex!, size, write, show

include("global.jl")
export DEFAULT_AA_MAPPING

include("objects.jl")
export DCAGraph, DCASample, eachsequence


include("IO.jl")
export writeparam

include("alignment_tools.jl")
export read_msa_num, read_msa

include("alignment_frequencies.jl")
export pairwise_frequencies, computeweights

include("modeltools.jl")
export switchgauge!, energy

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
