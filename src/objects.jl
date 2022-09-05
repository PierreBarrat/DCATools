"""
	mutable struct DCAGraph

Structure representing a Potts model.

## Fields
```
L::Int64
q::Int64
J::Array{Float64,2}
h::Array{Float64,1}
```
"""
Base.@kwdef mutable struct DCAGraph
    L::Int64
    q::Int64
    J::Array{Float64,2} = zeros(Float64, L*q, L*q)
    h::Array{Float64,1} = zeros(Float64, L*q)
end
"""
	DCAGraph(L,q)

Creates a `DCAGraph` of the required size. Parameters are initialized to zero.
"""
function DCAGraph(L,q)
	return DCAGraph(zeros(Float64, L*q, L*q), zeros(Float64, L*q), L, q)
end


Base.size(g::DCAGraph) = (g.q, g.L)

"""
	*(B, g::DCAGraph)

Multiply fields and couplings in `g` by scalar `B`. Useful to change temperature.
"""
function Base.:*(B, g::DCAGraph)
    return DCAGraph(B*g.J, B*g.h, g.L, g.q)
end

"""
	*(B, g::DCAGraph)

Multiply fields and couplings in `g` by scalar `B`. Useful to change temperature.
"""
function Base.:*(g::DCAGraph, B)
    return DCAGraph(B*g.J, B*g.h, g.L, g.q)
end

"""
"""
Base.getindex(g::DCAGraph, i, j, a, b) = g.J[(i .-1)*g.q .+ a, (j .-1)*g.q .+ b]
Base.getindex(g::DCAGraph, i, j, a, b::Colon) = g.J[(i .-1)*g.q .+ a, (j .-1)*g.q .+ (1:g.q)]
Base.getindex(g::DCAGraph, i, j, a::Colon, b) = g.J[(i .-1)*g.q .+ (1:g.q), (j .-1)*g.q .+ b]
Base.getindex(g::DCAGraph, i, j, a::Colon, b::Colon) = g.J[(i .-1)*g.q .+ (1:g.q), (j .-1)*g.q .+ (1:g.q)]

Base.getindex(g::DCAGraph, i, a) = g.h[(i .-1)*g.q .+ a]
Base.getindex(g::DCAGraph, i, a::Colon) = g.h[(i .-1)*g.q .+ (1:g.q)]
Base.getindex(g::DCAGraph, i::Colon, a) = g.h[(0:g.L-1)*g.q .+ a]
Base.getindex(g::DCAGraph, i::Colon, a::Colon) = g.h[:]

Base.setindex!(g::DCAGraph, val, i, j, a, b) = (g.J[(i .-1)*g.q .+ a, (j .-1)*g.q .+ b] = val)
Base.setindex!(g::DCAGraph, val, i, j, a::Colon, b) = (g.J[(i .-1)*g.q .+ (1:g.q), (j .-1)*g.q .+ b] .= val)
Base.setindex!(g::DCAGraph, val, i, j, a, b::Colon) = (g.J[(i .-1)*g.q .+ a, (j .-1)*g.q .+ (1:g.q)] .= val)
Base.setindex!(g::DCAGraph, val, i, j, a::Colon, b::Colon) = (g.J[(i .-1)*g.q .+ (1:g.q), (j .-1)*g.q .+ (1:g.q)] .= val)


Base.setindex!(g::DCAGraph, val, i, a) = (g.h[(i .-1)*g.q .+ a] = val)
Base.setindex!(g::DCAGraph, val, i, a::Colon) = (g.h[(i .-1)*g.q .+ (1:g.q)] .= val)
Base.setindex!(g::DCAGraph, val, i::Colon, a) = (g.h[(0:g.L-1)*g.q .+ a] .= val)
Base.setindex!(g::DCAGraph, val, i::Colon, a::Colon) = (g.h[:] .= val)
