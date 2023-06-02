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
    mapping::String = DEFAULT_MAPPING(q)
end

function DCAGraph(
	J::AbstractArray{<:Real, 4}, h::AbstractArray{<:Real, 2};
	q = size(J,3), mapping = DEFAULT_MAPPING(q),
)
	@assert size(J,3) == size(J,4) == size(h,2) "Incoherent sizes for J and h" size(J) size(h)
	@assert size(J,1) == size(J,2) == size(h,1) "Incoherent sizes for J and h" size(J) size(h)

	L = size(J,1)

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

	return DCAGraph(L, q, J_, h_, mapping)
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
	mapping = DEFAULT_AA_MAPPING,
)
	if init == :null
		DCAGraph(L, q, zeros(Float64, L*q, L*q), zeros(Float64, L*q), mapping)
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


"""
	mutable struct DCASample

```
    dat::Matrix{Int}
    q::Int = 21
    mapping::String = DCATools.DEFAULT_MAPPING(q)
    weights::Vector{Float64} = ones(size(dat,1))/size(dat,1) # sums to 1
```

Stores sequences or a sample of a DCA model. `dat` stores sequences/samples in *columns*:
`eachcol(X.dat)` will iterate over sequences.

**Important**: When built from a matrix, will *transpose* the input; if `size(dat) = (M, L)`,
`X=DCASample(dat)` will return an object with `size(X.dat) = (L, M)`. In other words, assumes
that the input matrix has sequences as rows.

## Methods

- `getindex(X::DCASample, i)` returns a matrix/vector `X.dat[:, i]`.
- `for s in X::DCASample` iterates over sequences.
- `eachsequence(X::DCASample)` returns an iterator over sequences (`Vector{Int}`).
- `eachsequence_weighted(X::DCASample)` returns an iterator over sequences and weights.
- `subsample(X::DCASample, i)` constructs the subsample defined by index `i`.

"""
Base.@kwdef mutable struct DCASample
	dat::Matrix{Int}
	q::Int = 21
	mapping::String = DEFAULT_MAPPING(q)
	weights::Vector{Float64} = ones(size(dat,1))/size(dat,1)
	function DCASample(dat, q, mapping, weights)
		@assert isempty(mapping) || q == length(mapping) "Inconsistent size for mapping $mapping and q=$q.
		Use `mapping=\"\"` if you do not care about the mapping."
		@assert length(weights) == size(dat, 1) "inconsistent number of weights"
		@assert all(map(>(0), weights)) "Weights cannot be negative"
		@assert isapprox(sum(weights), 1) "Weights must sum to 1"
		new(Matrix(dat'), q, mapping, weights)
	end
end

"""
    DCASample(Y; q=21, kwargs...)
    DCASample(Y, q; kwargs...)

Build a sample from matrix or vector `Y`. Assume that `Y` has sequences as rows.

Keyword arguments are the fields of `DCASample`.
"""
DCASample(Y::AbstractMatrix; kwargs...) = DCASample(dat=Y; kwargs...)
DCASample(Y::AbstractMatrix, q; kwargs...) = DCASample(Y; q, kwargs...)
DCASample(s::AbstractVector; kwargs...) = DCASample(s[:,:]'; kwargs...)
DCASample(s::AbstractVector, q; kwargs...) = DCASample(s; q, kwargs...)

Base.iterate(X::DCASample) = iterate(eachcol(X.dat))
Base.iterate(X::DCASample, state) = iterate(eachcol(X.dat), state)
Base.eltype(::Type{DCASample}) = AbstractVector{Int}
Base.size(X::DCASample) = size(X.dat)
Base.size(X::DCASample, i::Int) = size(X.dat, i)
Base.length(X::DCASample) = size(X, 2)

Base.getindex(X::DCASample, i) = X.dat[:, i]
Base.firstindex(::DCASample) = 1
Base.lastindex(X::DCASample) = length(X)
Base.view(X::DCASample, i) = view(X.dat, :, i)



function Base.show(io::IO, X::DCASample)
	M, L = size(X)
	print(io, "Alignment of $M sequences of length $L - ")
	show(io, X.dat')
end
function Base.show(io::IO, x::MIME"text/plain", X::DCASample)
	L, M = size(X)
	println(io, "Alignment of $M sequences of length $L")
	show(io, x, X.dat')
end

eachsequence(X::DCASample) = eachcol(X.dat)
eachsequence_weighted(X::DCASample) = zip(eachsequences(X), X.weights)

function subsample(X::DCASample, i::Int)
    dat = reshape(X[i], length(X[i]))
    w = [1]
    return DCASample(dat; q=X.q, mapping=X.mapping, weights = w)
end
function subsample(X::DCASample, idx)
    dat = X[idx]'
    w = X.weights[idx]
    w /= sum(w)
    return DCASample(dat, X.q; mapping=X.mapping, weights = w)
end
