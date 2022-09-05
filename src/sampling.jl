"""
	sample(
		graph::DCAGraph,
		M;
		init = rand(1:graph.q, graph.L),
		Twait = 1,
		burnin = 5*Twait,
		beta = 1.0,
		outfile = "",
		rng = Random.GLOBAL_RNG,
		verbose = false
	)

Sample `M` configurations from probability distribution defined by `graph`.

## Kwargs
- `init`: Initial configuration.
- `Twait`: Number of MCMC *sweeps* between two samples.
- `burnin`: Burnin before first sample is taken.
- `beta`: Inverse temperature.
- `outfile`: File path to save sample as a csv. Empty string to not save.
"""
function sample(
	graph::DCAGraph,
	M;
	init = rand(1:graph.q, graph.L),
	Twait = 1,
	burnin = 5*Twait,
	beta = 1.0,
	outfile = "",
	rng = Random.GLOBAL_RNG,
	verbose = false
)
	# Unnecessary allocation here?
    if beta != 1.
		graph_local = deepcopy(graph)
        graph_local.J .*= beta
        graph_local.h .*= beta
    else
    	graph_local = graph
    end


    # Initialisation
    verbose ? println("Initializing with $(burnin) burnin iterations... ") : print("")
    (q, L) = size(graph_local)
    conf = copy(init)
    # precomputed indices `(j-1)*q + conf[j]` - will be updated along `conf`
    jdx = convert.(UInt16, collect((j-1)*q+conf[j] for j in 1:L))

    # Burnin
    sample = zeros(Int, M, L)
    for t in 1:burnin
    	mcmc_sweep!(conf, jdx, graph_local; rng)
    end
    sample[1,:] = conf
    verbose ? println("done!") : print("")

    # Sampling
    verbose ? println("Sampling") : print("")
    for m = 2:M
        if verbose && mod(m, 500)==0
            @printf("It %d/%d           \r", m+1, M)
        end
        for t in 1:Twait
        	mcmc_sweep!(conf, jdx, graph_local; rng)
        end
        sample[m,:] .= conf
    end
    verbose ? println("Done") : print("")

    # Writing output
    if outfile!=""
    	writedlm(outfile, sample, ' ')
    end
    return sample
end


function mcmc_sweep!(conf, jdx, g; rng = Random.GLOBAL_RNG)
	q, L = size(g)
	@assert length(conf) == L "Configuration and graph have incompatible sizes\
		(respectively $(length(conf)) and $L"
	for rep in 1:L
		mcmc_step!(conf, jdx, g; rng)
	end
	return nothing
end

#=
Note: `jdx` contains precomputed indices for `conf`, of the form `(j-1)*q + conf[j]`
=#
function mcmc_step!(conf, jdx, g::DCAGraph; rng = Random.GLOBAL_RNG)
	E = 0.
	q, L = size(g)

	# Flip position
	i = rand(rng, 1:L)
	a = conf[i] # initial state
	b = mod(a - 1 + rand(rng, 1:(q-1)), q) + 1 # new state
	id_i_b = (i-1)*q+b # precomputed index for (i,b)
	# Delta E
	E = g.h[jdx[i]] - g.h[id_i_b]
	for j = 1:L
		if j != i
			E += g.J[jdx[i], jdx[j]] - g.J[id_i_b, jdx[j]]
		end
	end

	if E<=0. || exp(-E) > rand()
		conf[i] = b
		jdx[i] = (i-1)*q + b
	end
end




"""
	autocorr(sample::Array{Int64,2}, q::Int64)

Autocorrelation of each variable, corresponding to columns `sample`, as a function of time (*i.e.* lines of `sample). 
"""
function autocorr(sample::Array{Int64,2}, q::Int64)
	(M,L) = size(sample)
	navmin = M-25
	f1 = pairwise_frequencies(sample)[1]
	null = vec(sum(reshape(f1,q,L).^2, dims=1)')

	ac = zeros(Float64, M-navmin, L)
	for t in 1:(M-navmin)
		for m in 1:(M-t)
			for i in 1:L
				ac[t,i] += Float64(sample[m,i] == sample[m+t,i]) 
			end
		end
		ac[t,:] /= (M-t)
		ac[t,:] -= null
	end
	return ac
end

function estimatetau(g::DCAGraph ; itau = 1, M = 300, tol=0.05, nchains=5)

	S = [sample(g, M; Twait=itau) for n in 1:nchains]
	
	# Hamming distance vs time internally to a sample
	nav = 10
	h_intra = [DCATools.moving_average(hamming_v_time(S[n]), nav) / g.L for n in 1:nchains]
	
	# mean Hamming distance between samples 
	h_inter = maximum([mean(hamming_v_time(S1,S2)[M-100:end]) for S1 in S, S2 in S])/ g.L

	for tau in 1:length(h_intra[1])
		flag = true
		for n in 1:nchains
			if !isapprox(h_intra[n][tau], h_inter, atol=tol)
				flag = false
				break
			end
		end
		if flag
			return tau * itau, h_intra, h_inter
		end
	end
	return length(h_intra[1]) * itau, h_intra, h_inter
end


function hamming_v_time(sample)
	(M,L) = size(sample)
	navmin = 25
	h = zeros(Float64, M-navmin)
	for t in 1:(M-navmin)
		for m in 1:(M-t)
			h[t] += DCATools.hamming(sample[m,:], sample[m+t,:]) 
		end
		h[t] /= (M-t)
	end
	return h
end

function hamming_v_time(S1, S2)
	(M,L) = size(S1)
	if size(S2) != size(S1)
		@error "Samples of different sizes"
	end
	navmin = M-25
	h = Array{Float64,1}(undef, M)
	for m in 1:M
		h[m] = DCATools.hamming(S1[m,:], S2[m,:])
	end
	return h
end
