"""
    compute_mapping(s::AbstractString)

Return a `Dict{Int, Char}`: `Dict(i => c for (i,c) in enumerate(s))`.
"""
compute_mapping(s::AbstractString) = Dict(i => c for (i,c) in enumerate(s))

const _DEFAULT_MAPPING = "-ACDEFGHIKLMNPQRSTVWY"
const _MAPPING_Q2 = "01"
"""
	const DEFAULT_AA_MAPPING::Dict

`"-ACDEFGHIKLMNPQRSTVWY"`, as a `Dict{Int, Char}`.
"""
const DEFAULT_AA_MAPPING = compute_mapping(_DEFAULT_MAPPING)

"""
	default_mapping(q::Int)

- If `q==21`, return `DEFAULT_AA_MAPPING`.
- If `q<21`, return the restrition of `DEFAULT_AA_MAPPING` to the first `q` amino acids.
- If `q>21`, return a meaningless mapping `i => 'A'` for all `i`.
"""
function default_mapping(q::Int)
    @assert q > 0 "`q` must be strictly positive - got $q"

	return if q == 21
		DEFAULT_AA_MAPPING
    elseif q == 2
        compute_mapping(_MAPPING_Q2)
	elseif q < 21
		compute_mapping(_DEFAULT_MAPPING[1:q])
	else
        @warn "No clear defined mapping for q = $q > 21"
		Dict(i => 'A' for i in 1:q)
	end
end

aa_to_num(s; mapping = DEFAULT_AA_MAPPING) = aa_to_num(s, mapping)
function aa_to_num(s::AbstractString, mapping::AbstractDict)
    return map(x -> get(mapping, x, 1), collect(s)) # unknown chars automatically mapped to 1
end
"""
    num_to_aa(s::AbstractVector{Int}, mapping = DEFAULT_AA_MAPPING)
    num_to_aa(s; mapping = DEFAULT_AA_MAPPING)

Convert vector of integers to string using `mapping`.
"""
function num_to_aa(s::AbstractVector{<:Integer}, mapping::AbstractDict = DEFAULT_AA_MAPPING)
    return string(map(x -> mapping[x], s)...)
end
num_to_aa(s; mapping = DEFAULT_AA_MAPPING) = num_to_aa(s, mapping)
# function num_to_aa(S::DCASample)
#     return map(s -> num_to_aa(s; S.mapping), S)
# end
