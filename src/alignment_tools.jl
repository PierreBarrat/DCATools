"""
	read_msa_num(
        infile::AbstractString;
        index_style=1, header=false, mapping::Dict, compute_weights=false, theta, q,
    )

Read an MSA stored in `infile` in a numerical index_style. Return a `DCASample` object.
Kwargs `mapping` and `q` are passed to `DCASample`.

If `index_style=1`, amino acids should be mapped from 1 to `q`.
If `index_style=0`, they should be mapped from 0 to `q-1`.
`header=true` discards the first line of the input file.
"""
function read_msa_num(
	infile::AbstractString;
	index_style = 1,
	header = false,
	mapping = DEFAULT_AA_MAPPING,
	compute_weights = false,
	theta = 0.2,
	q = 0,
)
	Y = nothing
	try
		Y = if header
			readdlm(infile, Int, skipstart=1)
		else
			readdlm(infile, Int)
		end
	catch err
		println("inputoutput.jl - read_msa_num: readdlm failed to read $infile. The alignment may not be of the expected format.")
		error(err)
	end

	if index_style==0
		Y .+= 1
	elseif index_style==1
		if findmin(Y)[1] == 0
			error("inputoutput.jl - read_msa_num: file $infile contains a 0 value. `index_style` should be set to 0.")
		end
	end

	if q == 0
		q = maximum(Y)
	end

	weights = compute_weights ? computeweights(Y, theta) : ones(size(Y, 1))/size(Y,1)

	X = DCASample(Y, q, mapping, weights)
end

"""
	read_msa(fastafile::AbstractString; mapping = DEFAULT_AA_MAPPING)

Read a fasta file into a `DCASample` object.
"""
function read_msa(
	fastafile::AbstractString;
	mapping = DEFAULT_AA_MAPPING, compute_weights = false, theta = 0.2,
)
    @debug "Reading alignment with mapping $mapping"
    rev_mapping = Dict(c => i for (i,c) in mapping)
    Y = FASTA.Reader(open(fastafile, "r")) do fasta
        vcat(map(x -> aa_to_num(sequence(x), rev_mapping)', fasta)...)
    end
    names = FASTA.Reader(open(fastafile, "r")) do fasta
        map(description, fasta)
    end
	q = length(mapping)
	weights = compute_weights ? computeweights(Y, theta) : ones(size(Y, 1))/size(Y,1)

	return DCASample(Y, q, mapping, weights, names)
end



"""
    Base.write(file::AbstractString, S::DCASample; map=true)

Write `S` to `file`. If `map=false`, use a numerical format.
"""
function Base.write(file::AbstractString, S::DCASample; map=true, kwargs...)
	if map
		_write_fasta(file, S)
	else
		_write_num(file, S; kwargs...)
	end
end
function _write_fasta(file::AbstractString, S::DCASample)
    try
	    FASTAWriter(open(file, "w")) do io
	    	for (i,s) in enumerate(S)
	    		rec = FASTARecord("$i", DCATools.num_to_aa(s; mapping = S.mapping))
	    		write(io, rec)
	    	end
	    end
    catch err
        @warn "There was a problem when writing sequences to files;
        this could be due to an inadapted mapping, got $(S.mapping)."
    end
end
function _write_num(file::AbstractString, S::DCASample; header=false)
	open(file, "w") do io
		if header
			L = lenseq(S)
			write(io, "$L $(S.q)")
		end
		writedlm(io, S.dat', ' ')
	end
end




"""
    pw_hamming_distance(Y::DCASample; normalize=true, step=1)
    pw_hamming_distance(Y::DCASample, X::DCASample; normalize=true, step=1)

Pairwise hamming distance between sequences in `Y` and `X`. For `M` sequences, return a vector of
length `M(M-1)/2`. Only considers sequences at index `1:step:end` (useful for large
alignments).
"""
function pw_hamming_distance(Y::DCASample; normalize=true, step=1)
    L, M = size(Y)
    out = if step == 1
        Vector{Float64}(undef, Int(M*(M-1)/2))
    else
        N = sum(m1 -> length((m1+1):step:M), 1:step:M)
        Vector{Float64}(undef, N)
    end
    i = 1
    for m1 in 1:step:M, m2 in (m1+1):step:M
        out[i] = hamming(view(Y, m1), view(Y, m2))
        i += 1
    end

    return normalize ? out / L : out
end

function pw_hamming_distance(X::DCASample, Y::DCASample; normalize=true, step=1)
    L1, M1 = size(X)
    L2, M2 = size(Y)

    if L1 != L2
        throw(
            DimensionMismatch("samples must have the same sequence length: got $L1 != $L2")
        )
    end

    out = if step == 1
        Vector{Float64}(undef, M1*M2)
    else
        N = sum(m1 -> length(1:step:M1), 1:step:M2)
        Vector{Float64}(undef, N)
    end
    i = 1
    for m1 in 1:step:M1, m2 in 1:step:M2
        out[i] = hamming(view(X, m1), view(Y, m2))
        i += 1
    end

    return normalize ? out / L1 : out
end
