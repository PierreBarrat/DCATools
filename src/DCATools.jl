module DCATools

using DelimitedFiles
using Statistics
using LinearAlgebra
using FastaIO

import Base: *, getindex, ==

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

"""
"""
getindex(g::DCAgraph, i, j, a, b) = g.J[(i .-1)*g.q .+ a, (j .-1)*g.q .+ b]
getindex(g::DCAgraph, i, a) = g.h[(i .-1)*g.q .+ a]

function ==(g1::DCAgraph, g2::DCAgraph) 
    g1.J == g2.J && g1.h == g2.h && g1.L == g2.L && g1.q == g2.q
end



include("inputoutput.jl")
include("alignmenttools.jl")
include("modeltools.jl")
include("contactprediction.jl")
include("misc.jl")
include("MF.jl")

export DCAgraph, *

# Other sub-modules
include("MCMC.jl")
export doMCMC
include("MutLand.jl")
export SingleMut, Mutant, MutData, readmutdata, mapsinglemut!, findsinglemut, finddoublemut, computeepistasis!
include("BM.jl")
export bmlearn


end