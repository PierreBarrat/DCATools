"""
	read_msa_num(infile::AbstractString ; index_style=1, header=false, )

Read an MSA stored in `infile` in a numerical index_style.

If `index_style=1`, amino acids should be mapped from 1 to `q`.
If `index_style=0`, they should be mapped from 0 to `q-1`.
`header` argument allows for discarding the first line of `infile`.
"""
function read_msa_num(
	infile::AbstractString;
	index_style=1,
	header=false,
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
	mapping = DEFAULT_AA_MAPPING, compute_weights = false, theta = 0.2)
	mapdict = Dict(x=>findfirst(a->a==x, mapping) for x in mapping)

	fasta = FASTA.Reader(open(fastafile, "r"))
	Y = vcat(map(x -> map_aa_seq(sequence(x), mapdict)', fasta)...)
	q = length(mapping)
	weights = compute_weights ? computeweights(Y, theta) : ones(size(Y, 1))/size(Y,1)
	return DCASample(Y, q, mapping, weights)
end

map_aa_seq(s; mapping = DEFAULT_AA_MAPPING) = map_aa_seq(s, mapping)
function map_aa_seq(s::AbstractString, mapping::AbstractString)
	mapdict = Dict(x=>findfirst(a->a==x, mapping) for x in mapping)
	return map_aa_seq(s, mapdict)
end
function map_aa_seq(s::AbstractString, mapping::AbstractDict)
	return map(x -> get(mapping, x, 1), collect(s))
end
"""
	map_seq_to_aa(s::AbstractVector{Int}, mapping = DEFAULT_AA_MAPPING)
	map_seq_to_aa(s; mapping = DEFAULT_AA_MAPPING)
"""
function map_seq_to_aa(s::AbstractVector{Int}, mapping = DEFAULT_AA_MAPPING)
	return string(map(x -> mapping[x], s)...)
end
map_seq_to_aa(s; mapping = DEFAULT_AA_MAPPING) = map_seq_to_aa(s, mapping)

"""
    pdist(Y::Matrix{<:Integer})

Pairwise hamming distance between sequences in lines of `Y`.
"""
function pdist(Y::Matrix{<:Integer})
    M = size(Y,1)
    Yt = Y'
    out = zeros(Int64, M, M)
    for m1 in 1:M
        for m2 in (m1+1):M
            out[m2,m1] = hamming(Yt[:,m1], Yt[:,m2])
        end
    end
    return out + out'
function write(file::AbstractString, S::DCASample; map=true, kwargs...)
	if map
		_write_fasta(file, S)
	else
		_write_num(file, S; kwargs...)
	end
end
function _write_fasta(file::AbstractString, S::DCASample)
	FASTAWriter(open(file, "w")) do io
		for (i,s) in enumerate(S)
			rec = FASTARecord("$i", DCATools.map_seq_to_aa(s; mapping = S.mapping))
			write(io, rec)
		end
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
