## This could be struct instead of mutable struct?
"""
	mutable struct DCAGraph

Structure representing a Potts model.
Can be constructed with keywords, *e.g.* `DCAGraph(; L=2, q=2, J=..., h=...)`.

## Fields
```
L::Int
q::Int
J::Array{Float64,2}
h::Array{Float64,1}
```
"""
Base.@kwdef mutable struct DCAGraph
    L::Int
    q::Int
    J::Array{Float64,2} = zeros(Float64, L*q, L*q)
    h::Array{Float64,1} = zeros(Float64, L*q)
end

function DCAGraph(J::AbstractArray{<:Real, 4}, h::AbstractArray{<:Real, 2})
	@assert size(J,3) == size(J,4) == size(h,2) "Incoherent sizes for J and h" size(J) size(h)
	@assert size(J,1) == size(J,2) == size(h,1) "Incoherent sizes for J and h" size(J) size(h)

	L = size(J,1)
	q = size(J,3)

	J_ = zeros(L*q, L*q)
	for i in 1:L, j in i:L
		J_[(i-1)*q .+ (1:q), (j-1)*q .+ (1:q)] .= J[i,j,:,:]
		J_[(j-1)*q .+ (1:q), (i-1)*q .+ (1:q)] .= J[j,i,:,:]
	end
	@assert J_ == J_' "`J` matrix is not symetric"

	h_ = zeros(L*q)
	for i in 1:L
		h_[(i-1)*q .+ (1:q)] .= h[i,:]
	end

	return DCAGraph(L, q, J_, h_)
end

"""
	DCAGraph(
		L,q;
		init = :null, Jrand = N -> 1/L*randn(N,N), hrand = N -> 1/sqrt(L)*randn(N),
	)

Return a `DCAGraph` of the required size.
- `init == :null`: parameters are intialized to zero.
- `init == :rand`: parameters are randomly sampled using `Jrand` and `hrand` keywords.

## Random initialization

- `hrand` should be a function `N::Int -> h::Vector{Float64}` with `h` of size `N`.
- `Jrand` should be a function `N::Int -> J::Matrix{Float64}` with `J` of size `NxN`.

`Jrand` does not have to return a symetric matrix.
The output matrix is made symetric with zeroes on the diagonal blocks.

"""
function DCAGraph(
	L, q;
	init = :null, Jrand = N -> 1/L*randn(N,N), hrand = N -> 1/sqrt(L)*randn(N),
)
	if init == :null
		DCAGraph(L, q, zeros(Float64, L*q, L*q), zeros(Float64, L*q))
	elseif init == :rand
		return random_graph(L, q, Jrand, hrand)
	end
end

function random_graph(L, q, Jrand, hrand)
	h = hrand(L*q)
	J = Jrand(L*q)
	J .= (J .+ J')/2
	for i in 1:L
		J[(i-1)*q .+ (1:q), (i-1)*q .+ (1:q)] .= 0
	end
	return DCAGraph(; L, q, J, h)
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
