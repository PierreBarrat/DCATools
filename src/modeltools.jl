"""
    switchgauge!(g::DCAGraph; gauge=:sum0, col=g.q, wt)

Switch parameters in `g` to gauge `gauge`.
Implemented gauges:
- 0 sum: `:sum0`
- Lattice gas (state `col` has energy 0): `:LG`, `:lg`, `:latticegas`.
- Wild type (sequence `wt` has energy 0): `:wt`
"""
function switchgauge!(g::DCAGraph; gauge=:sum0, col=g.q, wt=Array{Int,1}(undef,0))
    if gauge==:sum0
        g.J,g.h = switchgaugesum0(g.J,g.h,g.L,g.q)
    elseif gauge==:LG || gauge==:lg ||  gauge==:latticegas
        g.J,g.h = switchgaugeLG(g.J,g.h,g.q,col=col)
    elseif gauge==:wt
        g.J, g.h = switchgaugeWT(g.J, g.h, wt, g.q)
    elseif gauge==:ML || gauge==:ml
        g.J, g.h = switchgaugeML(g.J, g.h, g.L, g.q)
    else
        error("modeltools.jl - switchgauge: unrecognized `gauge=$(gauge)` keyword.")
    end
    return g
end

"""
    switchgaugesum0(J, h, L, q)

Switch parameters to 0 sum gauge.
"""
function switchgaugesum0(J, h, L, q)
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
    switchgaugeLG(J, h, q ; col=q)

Switch parameters to lattice gas gauge. Keyword `col` indicates the state for which parameters should be 0. Defaults to `q`. 
"""
function switchgaugeLG(J, h, q; col=q)
    L = Int(size(J,1)/q)
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
    switchgaugeWT(J, h, wt::Array{Int,1})
"""
function switchgaugeWT(J, h, wt, q)
    L = Int(size(J,1)/q)
    ho = zeros(Float64, L*q)
    Jo = zeros(Float64, L*q, L*q)
    t = zeros(Float64,q,q)

    ho .= h
    for i in 1:L
        for j in 1:L
            t = J[(i-1)*q .+ (1:q),(j-1)*q .+ (1:q)]
            Jo[(i-1)*q .+ (1:q),(j-1)*q .+ (1:q)] = t - repeat(t[:,wt[j]],1,q) - repeat(t[[wt[i]],:],q,1) + repeat([t[wt[i],wt[j]]],q,q)
            ho[(i-1)*q .+ (1:q)] +=  vec(t[:,wt[j]]) - vec(repeat([t[wt[i],wt[j]]],1,q))
        end
        ho[(i-1)*q .+ (1:q)] -= vec(repeat([h[(i-1)*q+wt[i]]],1,q))
    end

    return (Jo,ho)    
end

"""
"""
function switchgaugeML(J, h, L, q)
    Jo, ho = switchgaugesum0(J,h,L,q)
    α = 1. /(q + L - 1)
    for i in 1:L
        for j in (i+1):L
            Jo[(i-1)*q .+ (1:q), (j-1)*q .+ (1:q)] += α*repeat(ho[(i-1)*q .+ (1:q)], 1, 21) + α*repeat(ho[(j-1)*q .+ (1:q)]', 21, 1)
            Jo[(j-1)*q .+ (1:q), (i-1)*q .+ (1:q)] .= Jo[(i-1)*q .+ (1:q), (j-1)*q .+ (1:q)]'
        end
        ho[(i-1)*q .+ (1:q)] *= α*q 
    end
    return Jo, ho
end

"""
    energy(g::DCAGraph, sample)

Compute energies of all configurations in `sample` with graph `g`.
"""
function energy(g::DCAGraph, S::AbstractVector)
	L = length(S)
    L != g.L && throw(ErrorException("Incorrect sequence length $L. Model length $(g.L)"))
	E = 0
	for i in 1:L
        for j in (i+1):L
            E -= g[j, i, S[j], S[i]]
        end
    	E -= g[i, S[i]]
    end
    return E
end
energy(g::DCAGraph, S::AbstractMatrix{Int}) = map(s -> energy(g,s), eachrow(S))
energy(g::DCAGraph, S::DCASample) = map(s -> energy(g, s), S)


"""
    profile_model(X::DCASample; pc = 1e-5, as_graph = false)

Infer a  `ProfileModel` from alignment `X`. If `as_graph`, return a `DCAGraph` with zero
couplings instead.
"""
function profile_model(X::DCASample; pc = 1e-5, as_graph = false)
    f1, _ = pairwise_frequencies(X)
    model = profile_model(f1, X.q)
    return as_graph ? DCAGraph(model; pc) : model
end

"""
    profile_model(f1::Array{Float64,1}; pc = 1e-5)

Infer a `ProfileModel` from frequencies `f1`.

Keywords:
- `pc`: Pseudocount ratio. Defaults to `1e-5`.
"""
function profile_model(w::Array{Float64,1}, q)
    L = Int(size(w,1)/q)
    return ProfileModel(; L, q, w)
end


"""
    pseudolikelihood(Y, g::DCAGraph; weights=ones(size(Y,1)))

Compute the pseudo-likelihood of configurations (*ie* sequences) in `Y` according to parameters in `g`. 
"""
function pseudolikelihood(
	Y, g::DCAGraph;
	weights=ones(Float64,size(Y,1))
)
    (M,L) = size(Y)
    Yi = Y'

    PL = 0.
    Z = 0.
    p = 0.
    E = 0.
    for m in 1:M
        for i in 1:L
            # Local partition function
            Z = 0.
            for a in 1:g.q
                E = 0.
                for j in 1:L
                    E += g.J[(j-1)*g.q+Yi[j,m], (i-1)*g.q+a] + g.h[(j-1)*g.q+Yi[j,m]]
                end
                Z += exp(E)
            end
            # Local likelihood of data
            E = 0.
            for j in 1:L
                E += g.J[(j-1)*g.q+Yi[j,m], (i-1)*g.q+Yi[i,m]] + g.h[(j-1)*g.q+Yi[j,m]]
            end
            PL += log(exp(E)/Z) * weights[m]
        end
    end
    PL = PL/sum(weights)
    return PL
end

