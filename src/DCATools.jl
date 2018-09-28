module DCATools

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

import Base: *
"""
	function *(B, g::DCAgraph)

Multiply fields and couplings in `g` by scalar `B`. Useful to change temperature. 
"""
function *(B, g::DCAgraph)
    return DCAgraph(B*g.J, B*g.h, g.L, g.q)
end

"""
	function *(B, g::DCAgraph)

Multiply fields and couplings in `g` by scalar `B`. Useful to change temperature. 
"""
function *(g::DCAgraph, B)
    return DCAgraph(B*g.J, B*g.h, g.L, g.q)
end


include("inputoutput.jl")
include("alignmenttools.jl")
include("modeltools.jl")
include("contactprediction.jl")

export DCAgraph


end