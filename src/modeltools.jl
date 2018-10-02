export switchgauge!, computeenergies, inferprofile

"""
    switchgauge!(g::DCAgraph ; gauge="0sum")

Switch parameters in `g` to gauge `gauge`. 
Implemented gauges:
1. 0 sum: "0sum"
2. Lattice gas: "LG", "lg", "latticegas". 
"""
function switchgauge!(g::DCAgraph ; gauge="0sum", col=g.q)
    if gauge=="0sum"
        g.J,g.h = switchgauge0sum(g.J,g.h,g.L,g.q)
    elseif gauge=="LG" || gauge=="lg" ||  gauge=="latticegas"
        g.J,g.h = switchgaugeLG(g.J,g.h,g.q,col=col)
    else
        error("modeltools.jl - switchgauge: unrecognized `gauge` keyword.")
    end
    return g
end

"""
    switchgauge0sum(J::Array{Float64,2}, h::Array{Float64,1}, L::Int64, q::Int64)

Switch parameters to 0 sum gauge.
"""
function switchgauge0sum(J::Array{Float64,2}, h::Array{Float64,1}, L::Int64, q::Int64)
    ho = zeros(Float64, L*q)
    Jo = zeros(Float64, L*q, L*q)
    t = zeros(Float64,q,q)
    
    for i in 1:L
        for j in (i+1):L
            t = J[(i-1)*q .+ (1:q), (j-1)*q .+ (1:q)]
            Jo[(i-1)*q .+ (1:q), (j-1)*q .+ (1:q)] = t - repeat(Statistics.mean(t,dims=1),q,1) - repeat(Statistics.mean(t,dims=2),1,q) .+ Statistics.mean(t)
        end
    end 
    Jo = Jo + Jo'

    for i in 1:L
        ho[(i-1)*q .+ (1:q)] = h[(i-1)*q .+ (1:q)] .- Statistics.mean(h[(i-1)*q .+ (1:q)])
        for j in 1:L
            t = J[(i-1)*q .+ (1:q), (j-1)*q .+ (1:q)]
            ho[(i-1)*q .+ (1:q)] = ho[(i-1)*q .+ (1:q)] + Statistics.mean(t,dims=2) .- Statistics.mean(t)
        end
    end
    return (Jo,ho)
end

"""
    switchgaugeLG(J::Array{Float64,2}, h::Array{Float64,1}, q::Int64 ; col=q)

Switch parameters to lattice gas gauge. Keyword `col` indicates the state for which parameters should be 0. Defaults to `q`. 
"""
function switchgaugeLG(J::Array{Float64,2}, h::Array{Float64,1}, q::Int64 ; col=q)
    L = Int64(size(J,1)/q)
    ho = zeros(Float64, L*q)
    Jo = zeros(Float64, L*q, L*q)
    t = zeros(Float64,q,q)

    ho .= h
    for i in 1:L
        for j in 1:L
            t = J[(i-1)*q .+ (1:q),(j-1)*q .+ (1:q)]
            Jo[(i-1)*q .+ (1:q),(j-1)*q .+ (1:q)] = t - repeat(t[:,q],1,q) - repeat(t[[q],:],q,1) + repeat([t[q,q]],q,q)
            ho[(i-1)*q .+ (1:q)] +=  vec(t[:,q]) - vec(repeat([t[q,q]],1,q))
        end
        ho[(i-1)*q .+ (1:q)] -= vec(repeat([h[(i-1)*q+q]],1,q))
    end

    return (Jo,ho)
end


"""
    computeenergies(g::DCAgraph, sample::Array{Int64,2})

Compute energies of all configurations in `sample` with graph `g`.
"""
function computeenergies(g::DCAgraph, sample::Array{Int64,2})

    (M,L) = size(sample)
    energies = zeros(Float64, M)
    for m = 1:M
        for i = 1:L
            for j = (i+1):L
                energies[m] -= g.J[(i-1)*g.q+sample[m,i], (j-1)*g.q+sample[m,j]]
            end
        energies[m] -= g.h[(i-1)*g.q+sample[m,i]]
        end
    end
    return energies
end

"""
    computeenergies(g::DCAgraph, sample::Array{Int64,1})

Compute energies of all configurations in `sample` with graph `g`.
"""
function computeenergies(g::DCAgraph, sample::Array{Int64,1})
    return computeenergies(g,sample[:])
end


"""
    inferprofile(Y::Array{Int64,2}; q=findmax(Y)[1], pc = 1e-5, weights=[], save::String="")

Infer profile model from alignment `Y`. 

Keywords:
- `q`: Default to maximum value in `Y`. 
-`pc` and `save`: See `inferprofile`. 
- `weights`: see `computefreqs`
- 
"""
function inferprofile(msa::Array{Int64,2}, q::Int64; pc = 1e-5, weights="", outfile::String="")
    f1 = computefreqs(msa, q=q, weights=weights)
    return inferprofile(f1, q, pc=pc, save=save)
end

"""
    inferprofile(f1::Array{Float64,1}, q::Int64; pc = 1e-5, save::String="")

Infer profile model from frequencies `f1`. 

Keywords:
-`pc`: Pseudocount ratio. Defaults to `1e-5`.
- save: File to save inferred profile. Format of save is `mat`, readable by `readdlm` or `readparam`. 
"""
function inferprofile(f1::Array{Float64,1}, q::Int64; pc = 1e-5, save::String="")
    L = Int64(size(f1,1)/q)
    h = log.((1-pc)*f1 .+ pc/q)
    for i in 1:L
        h[(i-1)*q+(1:q)] .-= mean(h[(i-1)*q+(1:q)])
    end
    if save!=""
        writedlm(outfile, [zeros(L*q,L*q) ; h'], " ")
    end
    return DCAgraph(zeros(L*q,L*q),h,L,q)
end

