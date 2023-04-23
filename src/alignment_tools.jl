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
end

"""
	read_msa_num(infile::AbstractString ; index_style=1, header=false)

Read an MSA stored in `infile` in a numerical index_style.

If `index_style=1`, amino acids should be mapped from 1 to `q`.
If `index_style=0`, they should be mapped from 0 to `q-1`.
`header` argument allows for discarding the first line of `infile`.
"""
function read_msa_num(infile::AbstractString ; index_style=1, header=false)
	Y = Array{Float64,2}(undef,0,0)
	try
		if header
			Y = readdlm(infile, Int, skipstart=1)
		else
			Y = readdlm(infile, Int)
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

	return Y
end

# """
#     convert_fasta(infasta::String, outfasta::String, mapping)

# Convert amino-acid characters in `infasta` to numbers, using `mapping.
# Write result to `outfasta`.
# """
# function convert_fasta(infasta::String, outfasta::String, mapping = DEFAULT_AA_MAPPING)
#     fasta = readfasta(infasta)
#     out = Array{Int64,2}(undef, length(fasta), length(fasta[1][2]))
#     for (i,(n,s)) in enumerate(fasta)
#         out[i,:] .= convert_seq(s, mapping)
#     end
#     writedlm(outfasta, out, ' ')
# end

# function convert_seq(s::String, mapping = DEFAULT_AA_MAPPING)
#     mapdict = Dict(x=>findfirst(a->a==x, mapping) for x in mapping)
#     σ = Array{Int64,1}(undef, length(s))
#     for (i,a) in enumerate(s)
#         σ[i] = get(mapdict, a, 1)
#     end
#     return σ
# end
