"""
    fitquality(f2_1::Array{Float64,2}, f1_1::Array{Float64,1}, f2_2::Array{Float64,2}, f1_2::Array{Float64,1}, q::Int64; withdiag=false)

Fitting quality of pairwise frequencies. Compare frequencies `(f1_1,f2_1)` (*e.g.* from natural sequences) to `(f1_2,f2_2)` (*e.g.* from a dca model). Output is 
1. Pearson correlation between connected correlations. (default: without diagonal elements)
2. Slope corresponding to the Pearson correlation. 
3. Frobenius norm of the difference between connected correlations.
4. Same as 1. for magnetizations.
5. Same as 3. for magnetizations.
"""
function fitquality(f2_1::Array{Float64,2}, f1_1::Array{Float64,1}, f2_2::Array{Float64,2}, f1_2::Array{Float64,1}, q::Int64; withdiag=false)
    L = Int64(size(f2_1,1)/q)
    C1 = f2_1 - f1_1*f1_1'
    C2 = f2_2 - f1_2*f1_2'
    # id = zeros(Bool,L*q,L*q)

    c1vec = zeros(Float64, Int64(L*(L-1)*q*q/2) )
    c2vec = zeros(Float64, Int64(L*(L-1)*q*q/2) )
    pos = 0
    for i = 1:L
        # id[(i-1)*q .+ (1:q), (i-1)*q .+ (1:q)] .= withdiag
        for j = (i+1):L
            # id[(j-1)*q .+ (1:q),(i-1)*q .+ (1:q)] .= true
            c1vec[pos .+ (1:q^2)] .= vec(C1[(i-1)*q .+ (1:q), (j-1)*q .+ (1:q)])
            c2vec[pos .+ (1:q^2)] .= vec(C2[(i-1)*q .+ (1:q), (j-1)*q .+ (1:q)])
            pos += q^2
        end
    end

    cc = Statistics.cor(c1vec,c2vec)
    cslope = linreg(c1vec,c2vec)[2]
    froc = LinearAlgebra.norm(c1vec - c2vec)
    cm = Statistics.cor(f1_1,f1_2)
    from = LinearAlgebra.norm(f1_1-f1_2)
    return (cc,cslope,froc,cm,from)
end

"""
    linreg(x,y)

Simple linear regression. No longer available in Base sadly. 
"""
function linreg(x::Array{Float64,1},y::Array{Float64,1})
    return reverse([x ones(Float64, length(x))]\y)
end


## Three points correlation
"""
    mutable struct corr3p

Three body correlation at positions `i, j, k` for states `a, b, c`, with value `cor`. 
"""
mutable struct corr3p
    i::Int64
    j::Int64
    k::Int64
    a::Int64
    b::Int64
    c::Int64
    cor::Float64
end

"""
    threepointscor(Y::Array{Int64,2}, w::Array{Float64,1}, q::Int64, threshold::Float64)

Base routine. Compute three body correlations between columns of alignment `Y`, with sequences being weighted by `w`. Keep only correlation values that are above `threshold`. Return corresponding list of triplets `(i,j,k),(a,b,c)`. 
"""
function threepointscor(Y::Array{Int64,2}, w::Array{Float64,1}, q::Int64, threshold::Float64)
    (M,L) = size(Y)
    Meff = sum(w)
    @time (f1,f2) = pairwise_frequencies(Y,w,q)

    c3p = Array{corr3p,1}(undef, 0)
    ff3 = zeros(Float64, q*q*q) # Will contain 3p frequencies for triplet (i,j,k)

    @time for i in 1:L
    print("$i/$L ...    \r")
        for j in (i+1):L
             for k in (j+1):L
                # Computing 3p freqs
                ff3.=0
                for m in 1:M
                    ff3[(Y[m,i]-1)*q*q + (Y[m,j]-1)*q + Y[m,k]] += w[m] 
                end
                ff3/=Meff

                # Computing correlation
                for a in 1:q
                    for b in 1:q
                        for c in 1:q
                            ff3[(a-1)*q*q + (b-1)*q + c] -= f2[(i-1)*q+a, (j-1)*q+b]*f1[(k-1)*q+c] - f2[(i-1)*q+a, (k-1)*q+c]*f1[(j-1)*q+b] - f2[(j-1)*q+b, (k-1)*q+c]*f1[(i-1)*q+a] + 2*f1[(i-1)*q+a]*f1[(j-1)*q+b]*f1[(k-1)*q+c]
                            # Storing only if above threshold
                            if abs(ff3[(a-1)*q*q + (b-1)*q + c]) >= threshold
                                push!(c3p,corr3p(i,j,k,a,b,c,ff3[(a-1)*q*q + (b-1)*q + c]))
                                # triplets = [triplets ; [i,j,k,a,b,c]]
                            end
                        end
                    end
                end

            end
        end
    end

    # List of used triplets
    triplets = Array{Int64,2}(undef, size(c3p,1),6);
    for t in 1:size(c3p,1)
        triplets[t,:] = [c3p[t].i c3p[t].j c3p[t].k c3p[t].a c3p[t].b c3p[t].c]
    end

    return c3p, triplets
end

"""
    threepointscor(Y::Array{Int64,2}, w::Array{Float64,1}, q::Int64, triplets::Array{Int64,2})

Base routine. Compute three body correlations between columns of alignment `Y`, with sequences being weighted by `w`. Consider only triplets (i,j,k,a,b,c) in `triplets`. 

*Note*: this function could easily be optimized.
"""
function threepointscor(Y::Array{Int64,2}, w::Array{Float64,1}, q::Int64, triplets::Array{Int64,2})
    (M,L) = size(Y)
    Meff = sum(w)
    (f1,f2) = pairwise_frequencies(Y,w,q)

    K = size(triplets,1)
    c3p = Array{corr3p,1}(undef,0)
    ff3 = 0.
    tick = Int64(round(K/100))
    for tr in 1:K
        if mod(tr,tick)==0
            print("$(Int64(round(tr/K*100)))% complete    \r")
        end
        i = triplets[tr,1]
        j = triplets[tr,2]
        k = triplets[tr,3]
        a = triplets[tr,4]
        b = triplets[tr,5]
        c = triplets[tr,6]
        t = 0.
        for m in 1:M
            if Y[m,i] == a && Y[m,j] == b && Y[m,k] == c
                ff3+=w[m]
            end
        end
        ff3/=Meff
        ff3 -= f2[(i-1)*q+a, (j-1)*q+b]*f1[(k-1)*q+c] - f2[(i-1)*q+a, (k-1)*q+c]*f1[(j-1)*q+b] - f2[(j-1)*q+b, (k-1)*q+c]*f1[(i-1)*q+a] + 2*f1[(i-1)*q+a]*f1[(j-1)*q+b]*f1[(k-1)*q+c]
        push!(c3p,corr3p(i,j,k,a,b,c,ff3))
    end
    return c3p
end

"""
    threepointscor(Y::Array{Int64,2}; q = findmax(Y)[1], threshold=0, triplets = Array{Int64,2}(undef,0,1),computew=false, weights=[], theta=0.2, saveweights="")

Compute three points frequencies for an array input. `Y` is an array of `Int64`. Return an array of `c3p` structures. 

Keywords: 
- `q`: default `findmax(Y)[1]`
- `threshold`: default `0`. Only points with correlation higher than `threshold` are returned.
- `triplets`: default `Array{Int64,2}(undef,0,1)`. List of position and amino acid triplets `(i,j,k,a,b,c)` for which the correlation is computed. 
- `weights`: default `[]`. If it is a `String`, phylogenetic weights are read from the corresponding file. If it is an `Array{Float64,1}`, they are used directly.
- `computew`: default `false`. If true, phylogenetic weights are computed, calling the appropriate `computeweights`. `weights` is then ignored. 
- `saveweights` and `theta`: see `computeweights`. 
"""
function threepointscor(Y::Array{Int64,2}; q = findmax(Y)[1], threshold=0, triplets = Array{Int64,2}(undef,0,1), computew=false, weights=[], theta=0.2, saveweights="")
    w = Array{Float64,1}(undef, 0)
    if computew
        # compute weights
        w = computeweights(Y, theta=theta, saveweights=saveweights)
        if weights!=[]
            # computew cancels weights
            @warn("misc.jl - threepointscor: both keywords `weights` and `computew` were declared. `weights` ignored.\n")
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
            warn("misc.jl - threepointscor: unrecognized format for keyword `weights`.")
        end
    end
    if size(w,1)!=size(Y,1)
        error("misc.jl - threepointscor: incorrect size for `weights`. Should be equal to the number of sequences in `Y`.")
    end
    if triplets != Array{Int64,2}(undef,0,1)
        if threshold !=0
            @warn("misc.jl - threepointscor: both keywords `threshold` and `triplets` were declared. `threshold` ignored.\n")
        end
        c3p = threepointscor(Y,w,q,triplets)
    else
        c3p, triplets = threepointscor(Y,w,q,Float64(threshold))
    end
    return c3p, triplets
end





"""
    projseq(seq, pc, q)

Project sequences in `seq` on column vectors in `pc`. 
"""
function projseq(seq::Array{Int64,2}, pc::Array{Float64,2}, q)
    (K,P) = size(pc)
    (M,L) = size(seq)
    if K!=L*q
        error("misc.jl - projseq: incorrect size for `pc`.\n")
    end
    proj = zeros(Float64, M,P)
    seq01 = zeros(Bool,L*q)
    for m in 1:M
        seq01 .= false
        for i in 1:L
            seq01[(i-1)*q+seq[m,i]] = true
        end
        proj[m,:] = seq01' * pc
    end

    return proj
end

"""
    aa2bin(x::Array{Int64,1} ;q=21)
    aa2bin(x::Array{Int64,2} ;q=21)
"""
function aa2bin(x::Array{Int64,1} ;q=21)
    out = zeros(Int64, length(x)*q)
    for (i,a) in enumerate(x)
        out[(i-1)*q+a]=1
    end
    return out
end

function aa2bin(x::Array{Int64,2} ;q=21)
    out = zeros(Int64, size(x,1), size(x,2)*q)
    for m in 1:size(x,1)
        out[m,:] .= aa2bin(x[m,:], q=q)
    end
    return out
end

"""
    hamming(X,Y)
"""
function hamming(X::AbstractVector, Y::AbstractVector)
    @assert length(X) == length(Y) "Cannot compute hamming distance for arrays of different size"
    S = 0
    for (x,y) in zip(X,Y)
        if x != y
        	S += 1
        end
    end
    return S
end

function moving_average(X, n::Int64)
    if length(X) < n
        @error "`X` shorter than `n`=$n"
    end
    out = Array{Float64, 1}(undef, length(X)-n)
    for i in 1:length(X)-n
        # Av from X[i] to X[i+n]
        out[i] = mean(X[i:i+n])
    end
    return out
end

