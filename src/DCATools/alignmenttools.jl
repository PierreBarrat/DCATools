export computefreqs, computeweights, pdist, convert_fasta

"""
    computefreqs(Y::Array{Int64,2}, w::Array{Float64,1}, q::Int64)

Base routine for computing frequencies in sample `Y`. `w` is an array containing the weights. 
"""
function computefreqs(Y::Array{Int64,2}, w::Array{Float64,1}, q::Int64)
    if typeof(w)==Array{Float64,2}
        @warn("alignmenttools.jl - computefreqs: `w` is of dimension 2. Applying `vec`.")
        w = vec(w)
    end

    (M,L) = size(Y)
    f2 = zeros(Float64, L*q, L*q)
    f1 = zeros(Float64,L*q) 

    if size(w)[1]!=M
        error("alignmenttools.jl - computefreqs: incorrect number of weights\n")
    end 

    Meff = sum(w)
    for m in 1:M
        for j in 1:L
            f1[(j-1)*q+Y[m,j]] += w[m]
            f2[(j-1)*q+Y[m,j], (j-1)*q+Y[m,j]] += w[m]
            for i in (j+1):L
                f2[(i-1)*q+Y[m,i], (j-1)*q+Y[m,j]] += w[m]
                f2[(j-1)*q+Y[m,j], (i-1)*q+Y[m,i]] += w[m]
            end
        end
    end

    f2 = f2./Meff
    f1 = f1./Meff

    return (f1,f2)
end


"""
    computefreqs(Y::Array{Int64,2}; q = findmax(Y)[1], computew=false, weights=[], theta=0.2, saveweights="", pc=0.)

Compute pairwise frequencies for an array input. `Y` is an array of `Int64`. Return frequencies `f1` and `f2` and weights.

Keywords: 
- `q`: default `findmax(Y)[1]`
- `weights`: default `[]`. If it is a `String`, phylogenetic weights are read from the corresponding file. If it is an `Array{Float64,1}`, they are used directly.
- `computew`: default `false`. If true, phylogenetic weights are computed, calling the appropriate `computeweights`. `weights` is then ignored. 
- `saveweights` and `theta`: see `computeweights`. 
- `pc`: default 0. Pseudocount ratio. 
"""
function computefreqs(Y::Array{Int64,2}; q = findmax(Y)[1], computew=false, weights=[], theta=0.2, saveweights="", pc=0.)
    w = Array{Float64,1}(undef, 0)
    if computew
        # compute weights
        w = computeweights(Y, theta=theta, saveweights=saveweights)
        if weights!=[]
            # computew cancels weights
            @warn("alignmenttools.jl - computefreqs: both keywords `weights` and `computew` were declared. `weights` ignored.\n")
        end
    else
        if weights == []
            # no weights used
            w = ones(Float64, size(Y,1))
        elseif typeof(weights) == String
            # read them from file
            w = vec(readdlm(weights, Float64))
        elseif typeof(weights) == Array{Float64,1}
            # read them from the array
            w = weights
        else
            # not recognized, no weights used
            warn("alignmenttools.jl - computefreqs: unrecognized format for keyword `weights`.")
        end
    end
    if size(w,1)!=size(Y,1)
        error("alignmenttools.jl - computefreqs: incorrect size for `weights`. Should be equal to the number of sequences in `Y`.")
    end

    (f1, f2) = computefreqs(Y,w,q)
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

"""
    computefreqs(msa::String; q = 0, computew=false, weights=[], theta=0.2, saveweights="", header=false, format=1, pc=0.)

Compute pairwise frequencies for a file input. `msa` is a file containing the alignment in numerical format. Return frequencies `f1` and `f2` and weights.

Keywords: 
- `q`: default `findmax(Y)[1]`
- `weights`: default `[]`. If it is a `String`, phylogenetic weights are read from the corresponding file. If it is an `Array{Float64,1}`, they are used directly.
- `computew`: default `false`. If true, phylogenetic weights are computed, calling the appropriate `computeweights`. `weights` is then ignored. 
- `saveweights` and `theta`: see `computeweights`. 
- `header` and `format`: see `readmsanum`.
- `pc`: default 0. Pseudocount ratio.
"""
function computefreqs(msa::String; q = 0, computew=false, weights=[], theta=0.2, saveweights="", header=false, format=1, pc=0.)
    Y = readmsanum(msa, header=header, format=format)
    if q==0
        q = findmax(Y)[1]
    end
    return computefreqs(Y, q=q, computew=computew, weights=weights, theta=theta, saveweights=saveweights, pc=pc)
end



"""
    computeweights(msa::String; theta::Float64 = 0.2, saveweights::String = "", format=1, header=false)

Compute weights for file input. Compute weights of each sequence in file `msa` using reweighting threshold `theta` (default 0.2). 

Keywords: 
- `theta`: threshold of similarity under which sequences are weighted down. Default `0.2`. 
- `saveweights`: weights are saved there if non empty. Default `""`
- `format` and `header`: used to read `msa`. See `readmsanum`. 
"""
function computeweights(msa::String; theta::Float64 = 0.2, saveweights::String = "", format=1, header=false)
    Y = readmsanum(msa, format=format, header=header)
    w = computeweights(Y, theta=theta, saveweights=saveweights)
    return w
end

"""
    computeweights(Y::Array{Int64,2} ; theta = 0.2, saveweights="")  

Compute weights for array input. Compute weights of each sequence in file `msa` using reweighting threshold `theta` (default 0.2). 

Keywords: 
- `theta`: threshold of similarity under which sequences are weighted down. Default `0.2`. 
- `saveweights`: weights are saved there if non empty. Default `""`
"""
function computeweights(Y::Array{Int64,2} ; theta = 0.2, saveweights="") 
    w = computeweights(Y, theta)
    if saveweights!=""
        writedlm(saveweights,w," ")
    end
    return w
end

"""
    computeweights(Y::Array{Int64,2},theta::Float64)

Basic routine. Compute weights of sequences in alignment `Y`, using threshold `theta`. 
"""
function computeweights(Y::Array{Int64,2}, theta::Float64)
    Y = Y';
    M::Int64 = size(Y,2)
    L::Int64 = size(Y,1)
    h::Float64 = L * (1-theta)
    weights = ones(Float64, M);
    d::Int64=0
    for m = 1:M
        for l = (m+1):M
            d = 0
            for i = 1:L
                d += Y[i,m]==Y[i,l];
            end
            if d > h
                weights[m]+=1
                weights[l]+=1
            end
        end
    end
    return 1 ./weights
end

"""
    pdist(Y::Array{Int64,2})

Pairwise hamming distance between sequences in lines of `Y`. 
"""
function pdist(Y::Array{Int64,2})
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
    convert_fasta(infasta::String, outfasta::String)

Convert amino-acid characters in `infasta` to numbers, using the mapping `-ACDEFGHIKLMNPQRSTVWY`. 
Write result to `outfasta`.  
"""
function convert_fasta(infasta::String, outfasta::String)
    mapping = "-ACDEFGHIKLMNPQRSTVWY"
    mapdict = Dict(x=>findfirst(a->a==x, mapping) for x in mapping)
    fasta = readfasta(infasta)
    out = Dict()
    for (i,(n,s)) in enumerate(fasta)
        out[n] = ""
        for a in s
            out[n] *= "$(get(mapdict, a, "-")) "
        end
        out[n] = out[n][1:end-1]
    end
    out = map(x->parse.(Int64,x), split.(collect(values(out)), " "))
    out_ = zeros(Int64, length(fasta), length(out[1]))
    for (i,s) in enumerate(out)
        out_[i,:] .= s
    end
    # writefasta(outfasta, out)
    writedlm(outfasta, out_, ' ')
end