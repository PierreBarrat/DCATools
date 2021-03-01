module MCMC
	
using DCATools
using DelimitedFiles
using Printf
using Statistics
using Distributed
using Random


export doMCMC, samplefromgraph!, estimatetau, autocorr

"""
	doMCMC(graphfile::String,outfile::String, format::String, q::Int64,M::Int64,t ; T= 10000, beta = 1.0, verbose = true, conf_init=rand(1:q, L))

Takes file as input. Untested yet. 
"""
function doMCMC(graphfile::String,outfile::String, format::String, q::Int64,M::Int64,t ; T= 10000, beta = 1.0, verbose = true, conf_init=rand(1:q, L))
    g = readparam(graphfile, format=format, q=q)
    sample = doMCMC(g, M, t, outfile=outfile, T=T, beta=beta, verbose=verbose, conf_init=conf_init)
    return(sample)
end

"""
"""
function doMCMC_par(graph::DCAgraph, M::Int64, tau, nprocs::Int64; T= 50*tau, beta = 1.0, verbose = false)
	if length(workers()) < nprocs
		@error "Number of workers ($(length(workers()))) smaller than `nprocs=$nprocs`."
	end
	Mloc = cld(M, nprocs)
	if verbose
		println("Starting $nprocs MCMC chains of $Mloc samples each.")
	end
	sample = @distributed (vcat) for n in 1:nprocs
		doMCMC(graph, Mloc, tau, T = T, beta = beta, outfile="", verbose = false, conf_init = rand(1:graph.q, graph.L), nprocs = 1)
	end
	return sample
end

"""
	doMCMC(graph::DCAgraph, M::Int64, tau ; outfile="", T= 100000, beta = 1.0, verbose = false,conf_init=rand(1:graph.q, graph.L), nprocs = 1)

Sample `M` configurations from probability distribution defined by `graph`. Number of MCMC sweeps between configurations is `tau`. 

Keyword parameters: 
- outfile: Default "". Location to save sample. 
- T: Default 50*tau. Initial number of iterations. 
- beta: Default 1.0. Inverse temperature
- conf_init: Default random. Initial configuration. 
- nprocs: Default 1. Number of parallel MCMC chains.

## Note
Single swaps can be made by choosing `τ=1/L`
"""
function doMCMC(graph::DCAgraph, M::Int64, tau ; outfile="", T= 50*tau, beta = 1.0, verbose = false,conf_init=rand(1:graph.q, graph.L), nprocs = 1)
	# Handling parallel case first
	if nprocs > 1
		if verbose
			println("Parallel MCMC sampling with $nprocs workers.")
		end
		sample = doMCMC_par(graph, M, tau, nprocs, T = T, beta = beta, verbose=verbose)
    	if outfile!=""
    		writedlm(outfile,sample,' ')
    	end
		return sample
	end

	# Unnecessary allocation here?  
    if beta != 1.
		graph_local = deepcopy(graph)
        graph_local.J *= beta
        graph_local.h *= beta
    else
    	graph_local = graph
    end

    sample = zeros(Int64, M,graph_local.L)

    # Initialisation
    verbose ? println("Initializing with ",T," iterations... ") : print("")
    conf = rand(1:graph_local.q, graph_local.L)
    jdx = convert.(UInt16, collect((j-1)*graph_local.q+conf[j] for j in 1:graph_local.L))

    samplefromgraph!(conf, jdx, graph_local, T)
    sample[1,:] = conf
    verbose ? println("done!") : print("")

    # Sampling
    verbose ? println("Sampling") : print("")
    for m = 1:M-1
        if verbose && mod(m+1,500)==0
            @printf("It %d/%d           \r",m+1,M)
        end
        samplefromgraph!(conf, jdx, graph_local, tau)
        sample[m+1,:] .= conf
    end
    verbose ? println("Done") : print("")

    # Writing output
    if outfile!=""
    	writedlm(outfile,sample,' ')
    end
    return (sample)
end

"""
	samplefromgraph!(g::DCAgraph, conf_init::Array{Int64,1}, conf_end::Array{Int64,1}, tau::Int64)

Sample for `tau` sweeps from probability defined by `g`, starting with configuration `conf_init` and storing final configuration in `conf_end`. 
"""
function samplefromgraph!(conf::Array{Int64,1}, jdx::Array{UInt16,1}, g::DCAgraph, tau::Int64)
	E = 0.
	q = g.q
	L = g.L

	@fastmath @inbounds for t in 1:tau
		for i in 1:L
			a = conf[i]
			b = rand(1:g.q)
			while b==a
				b = rand(1:q)
			end

        	# id_i_a = (i-1)*q+a
        	id_i_b = (i-1)*q+b
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
	end
	nothing
end


"""
	samplefromgraph!(g::DCAgraph, conf_init::Array{Int64,1}, conf_end::Array{Int64,1}, tau::Float64)

Sample for `round(tau*L)` swaps from probability defined by `g`, starting with configuration `conf_init` and storing final configuration in `conf_end`. 

## Note
This version is used if `τ` is not an integer. It allows doing single swaps for `τ=1/L`  
"""
function samplefromgraph!(g::DCAgraph, conf_init::Array{Int64,1}, conf_end::Array{Int64,1}, tau::Float64)
	rng = MersenneTwister(rand(1:100000))
	E = 0.
	q = g.q
	L = g.L

	irng = Random.Sampler(rng, Set(1:L))
	brng = Random.Sampler(rng, Set(1:q))
	Erng = Random.Sampler(rng, Float64)
	copyto!(conf_end, conf_init)

	@fastmath @inbounds for t in 1:Int64(round(tau*L))
		i = rand(rng, irng)
		a = conf_end[i]
		b = rand(rng, brng)
		while b==a
			b = rand(rng, brng)
		end
        id_i_a = (i-1)*q+a
        id_i_b = (i-1)*q+b
        E = g.h[id_i_a] - g.h[id_i_b]
        for j = 1:L
        	if j != i
        		E += g.J[id_i_a, (j-1)*q+conf_end[j]] - g.J[id_i_b, (j-1)*q+conf_end[j]]
       		end
       	end
        if E<=0. || exp(-E) > rand(rng, Erng)
        	conf_end[i] = b
        end
    end
	nothing
end

"""
	estimatetau(g::DCAgraph)

Attempt to estimate reasonable number of iterations between samples. Based on autocorrelation. 
- Conservative: if the autocorrelation of the most autocorrelated spin is smaller than 1/e, eq. is reached. 
- Fast: if the average absolute autocorrelation etc... 
"""
# function estimatetau(g::DCAgraph ; itau = 1, M = 1000, threshold = 1/2.7, mode = "conservative", nprocs=1)

# 	t = doMCMC(g, M, itau, T=10000, nprocs=nprocs)
# 	ac = autocorr(t, g.q)
# 	out = Int64(0)

# 	if mode == "conservative"
# 		score = vec(findmax(ac, dims=2)[1])
# 	elseif mode == "fast"
# 		score = vec(mean(abs.(ac),dims=2))
# 	end
# 	if typeof(findnext(x->x<threshold, score, 1)) != Nothing
# 		out =  Int64(findnext(x->x<threshold, score, 1)) * itau
# 	else
# 		out =  size(ac,1) * itau
# 	end

# 	return Int64(out)
# end

"""
	autocorr(sample::Array{Int64,2}, q::Int64)

Autocorrelation of each variable, corresponding to columns `sample`, as a function of time (*i.e.* lines of `sample). 
"""
function autocorr(sample::Array{Int64,2}, q::Int64)
	(M,L) = size(sample)
	navmin = M-25
	f1 = computefreqs(sample)[1]
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

function estimatetau(g::DCAgraph ; itau = 1, M = 300, tol=0.05, nchains=5)

	S = [doMCMC(g, M, itau, T=100) for n in 1:nchains]
	
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



end # Module