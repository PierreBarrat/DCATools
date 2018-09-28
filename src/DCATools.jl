using DelimitedFiles
using Statistics

"""
	struct DCAgraph
"""
mutable struct DCAgraph
    J::Array{Float64,2}
    h::Array{Float64,1}
    L::Int64
    q::Int64
end

include("inputoutput.jl")
include("alignmenttools.jl")
include("modeltools.jl")
include("contactprediction.jl")