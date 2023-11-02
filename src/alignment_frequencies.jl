"""
	profile_frequency(X::DCASample)

Single site frequencies of `X`.
"""
function profile_frequency(X::DCASample)
	q = X.q
	(L, M) = size(X)

	f1 = zeros(Float64, L*q)
	for (m, s) in enumerate(X), (i, a) in enumerate(s)
		f1[(i-1)*q + a] += X.weights[m]
	end

	return f1
end
function consensus(X::DCASample)
    f = profile_frequency(X)
    L = size(X, 1)
    x = map(1:L) do i
        argmax(f[(i-1)*X.q .+ (1:X.q)])
    end
    return DCASample(x; mapping=X.mapping, names = ["consensus"])
end

"""
    pairwise_frequencies(Y::Matrix{<:Integer}, w::Array{Float64,1}, q::Integer)

Base routine for computing frequencies in sample `Y`. `w` is an array containing the weights.
"""
function pairwise_frequencies(Y, w::Array{Float64,1}, q)
    if typeof(w)==Array{Float64,2}
        @warn("`w` is of dimension 2. Applying `vec`.")
        w = vec(w)
    end

    (L, M) = size(Y)
    f2 = zeros(Float64, L*q, L*q)
    f1 = zeros(Float64, L*q)

    if length(w) != M
        error("incorrect number of weights\n")
    end
    Meff = sum(w)

    # off diagonal f2
    for m in 1:M, i in 1:L, j in (i+1):L
        f2[(i-1)*q + Y[i, m], (j-1)*q + Y[j, m]] += w[m]
	end
	f2 .= f2 .+ f2'

	# f1 and diagonal f2
    for m in 1:M, i in 1:L
    	f1[(i-1)*q + Y[i, m]] += w[m]
        f2[(i-1)*q + Y[i, m], (i-1)*q + Y[i, m]] += w[m]
    end
    f2 = f2 ./ Meff
    f1 = f1 ./ Meff

    return (f1, f2)
end


"""
    pairwise_frequencies(
    	Y
    	q = findmax(Y)[1], computew=false, weights=[], theta=0.2, saveweights="", pc=0.
    )

Compute pairwise frequencies for an array input. Columns of `Y` represent sequences.
Return frequencies `f1` and `f2` and weights.

## Kwargs:
- `weights`: default `[]`. If it is a `String`, phylogenetic weights are read from the \
	corresponding file. If it is an `Array{Float64,1}`, they are used directly.
- `computew`: default `false`. If true, phylogenetic weights are computed and \
	`weights` is ignored.
- `saveweights` and `theta`: see `computeweights`.
- `pc`: default 0. Pseudocount ratio.
"""
function pairwise_frequencies(
	Y::AbstractMatrix{<:Integer};
	q = findmax(Y)[1], computew=false, weights=[], theta=0.2, saveweights="", pc=0.
)
    w = Vector{Float64}(undef, 0)
    if computew
        # compute weights
        w = computeweights(Y, theta=theta, saveweights=saveweights)
        if weights!=[]
            # computew cancels weights
            @warn("both keywords `weights` and `computew` declared. `weights` ignored.\n")
        end
    else
        if weights == []
            # no weights used
            w = ones(Float64, size(Y,2))
        elseif typeof(weights) == String
            # read them from file
            w = vec(readdlm(weights, Float64))
        elseif typeof(weights) == Array{Float64,1}
            # read them from the array
            w = weights
        elseif typeof(weights) == Array{Float64,2}
            w = vec(weights)
        else
            # not recognized, no weights used
            @warn "unrecognized format for keyword `weights`"
        end
    end
    if length(w) != size(Y,2)
        error("incorrect size for `weights`: \
        	$(length(weights)) weights vs $(size(Y,2)) sequences")
    end

    (f1, f2) = pairwise_frequencies(Y, w, q)
    if pc != 0.
        f1 = (1-pc)*f1 .+ pc/q
        f2 = (1-pc)*f2 .+ pc/q/q
        for i in 1:size(Y,2)
            f2[(i-1)*q .+ (1:q), (i-1)*q .+ (1:q)] .= 0
            for a in 1:q
                f2[(i-1)*q+a, (i-1)*q+a] = f1[(i-1)*q+a]
            end
        end
    end
    return (f1,f2,w)
end

function pairwise_frequencies(Y::DCASample; computew=false, kwargs...)
	if computew
		f1, f2, w = pairwise_frequencies(Y.dat; q=Y.q, computew=true, kwargs...)
	else
		f1, f2, w = pairwise_frequencies(Y.dat; q=Y.q, weights=Y.weights, computew, kwargs...)
	end
	return f1, f2, w
end

"""
    pairwise_frequencies(
    	msa::String;
    	q = 0,
    	computew=false,
    	weights=[],
    	theta=0.2,
    	pc=0.,
    	msa_format=:numerical,
    	header = false,
    	saveweights="",
    )

Compute pairwise frequencies for the file input `msa`.
 Return frequencies `f1` and `f2` and weights.

## Kwargs
- `q`: default `findmax(Y)[1]`
- `weights`: default `[]`.
   If it is a `String`, phylogenetic weights are read from the corresponding file.
   If it is an `Array{Float64,1}`, they are used directly.
- `computew`: default `false`.
   If true, phylogenetic weights are computed, calling the appropriate `computeweights`.
   `weights` is then ignored.
- `saveweights` and `theta`: see `computeweights`.
- `header` and `msa_format`: see `read_msa_num`.
- `pc`: default 0. Pseudocount ratio.
"""
function pairwise_frequencies(
	msa::AbstractString;
	q = 0,
	pc = 0.,
	theta = 0.2,
	weights = [],
	computew = false,
	saveweights = "",
	msa_format = :numerical,
	index_style = 1,
	header = false,
	mapping = DEFAULT_AA_MAPPING,
)
	Y = if msa_format == :numerical
    	read_msa_num(msa; header, index_style, mapping)
	elseif msa_format == :fasta
		read_msa(msa; mapping)
	end

    return pairwise_frequencies(Y; Y.q, computew, weights, theta, saveweights, pc)
end



"""
    computeweights(
    	msa::AbstractString;
    	theta = 0.2, saveweights = "", msa_format=:numerical, header=false
    )

Compute weights for file input `msa` using threshold `theta`.

## Kwargs
- `theta`: threshold of similarity under which sequences are weighted down. Default `0.2`.
- `saveweights`: weights are saved there if non empty. Default `""`
- `msa_format` (`:fasta` or `:numerical`) and `header`: used to read `msa`. See `read_msa_num`.
"""
function computeweights(
	msa::AbstractString;
	theta = 0.2, saveweights = "", msa_format=:numerical, header=false
)
    Y = read_msa_num(msa, msa_format, header)
    w = computeweights(Y, theta, saveweights)
    return w
end

function computeweights!(Y::DCASample; theta=0.2)
	Y.weights = computeweights(Y.dat, theta)
	return Y.weights
end

"""
    computeweights(Y; theta = 0.2, saveweights="")

Compute phylogenetic weights for input `Y`.

## Kwargs
- `theta`: threshold of similarity under which sequences are weighted down. Default `0.2`.
- `saveweights`: weights are saved there if non empty. Default `""`
"""
function computeweights(Y::AbstractMatrix{<:Integer}; theta = 0.2, saveweights="")
    w = computeweights(Y, theta)
    if saveweights!=""
        writedlm(saveweights,w," ")
    end
    return w
end
computeweights(Y::DCASample; kwargs...) = computeweights(Y.dat; kwargs...)

"""
    computeweights(Y::Matrix{<:Integer}, theta::Float64)

Basic routine. Compute weights of sequences in alignment `Y`, using threshold `theta`.
"""
function computeweights(Y::AbstractMatrix{<:Integer}, theta::Float64)
    L, M = size(Y)
    h = L * (1-theta)
    weights = ones(Float64, M);
    d = 0
    for m = 1:M
        for l = (m+1):M
            d = 0
            for i = 1:L
            	if Y[i,m]==Y[i,l]
                	d += 1
                end
            end
            if d > h
                weights[m]+=1
                weights[l]+=1
            end
        end
    end
    weights = 1 ./ weights
    return weights / sum(weights)
end
