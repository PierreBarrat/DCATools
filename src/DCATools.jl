module DCATools

using DelimitedFiles
using Statistics
using LinearAlgebra
using FastaIO

"""
	DCAgraph
"""
mutable struct DCAgraph
    J::Array{Float64,2}
    h::Array{Float64,1}
    L::Int64
    q::Int64
end
function DCAgraph(L,q)
	return DCAgraph(zeros(Float64, L*q, L*q), zeros(Float64, L*q), L, q)
end

import Base: *

"""
	*(B, g::DCAgraph)

Multiply fields and couplings in `g` by scalar `B`. Useful to change temperature. 
"""
function *(B, g::DCAgraph)
    return DCAgraph(B*g.J, B*g.h, g.L, g.q)
end

"""
	*(B, g::DCAgraph)

Multiply fields and couplings in `g` by scalar `B`. Useful to change temperature. 
"""
function *(g::DCAgraph, B)
    return DCAgraph(B*g.J, B*g.h, g.L, g.q)
end


include("DCATools/inputoutput.jl")
include("DCATools/alignmenttools.jl")
include("DCATools/modeltools.jl")
include("DCATools/contactprediction.jl")
include("DCATools/misc.jl")

export DCAgraph, *


end