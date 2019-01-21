module DCAMCMC
	
using DCATools
using DelimitedFiles
using Printf
using Statistics
using Distributed
using Random
using Profile

export doMCMC, samplefromgraph!, estimatetau, autocorr

"""
	doMCMC(graphfile::String,outfile::String, format::String, q::Int64,M::Int64,t::Int64 ; T= 10000, beta = 1.0, verbose = true, conf_init=rand(1:q, L))

Takes file as input. Untested yet. 
"""
function doMCMC(graphfile::String,outfile::String, format::String, q::Int64,M::Int64,t::Int64 ; T= 10000, beta = 1.0, verbose = true, conf_init=rand(1:q, L))
    g = readparam(graphfile, format=format, q=q)
    sample = doMCMC(g, M, t, outfile=outfile, T=T, beta=beta, verbose=verbose, conf_init=conf_init)
    return(sample)
end

"""
"""
function doMCMC_par(graph::DCAgraph, M::Int64, tau::Int64, nprocs::Int64; T= 50*tau, beta = 1.0, verbose = false)
	if length(workers()) != nprocs
		@warn "Number of workers ($(length(workers()))) different from `nprocs=$nprocs` value."
	end
	Mloc = cld(M, nprocs)
	if Mloc != M/nprocs
		@warn "Sample size `M=$M` not divisible by `nproc=$nprocs`"
	end

	if verbose
		println("Starting $nprocs MCMC chains of $Mloc samples each.")
	end
	sample = @distributed (vcat) for n in 1:nprocs
		doMCMC(graph, Mloc, tau, T = T, beta = beta, outfile="", verbose = false, conf_init = rand(1:graph.q, graph.L), nprocs = 1)
	end
	return sample
end

"""
	doMCMC(graph::DCAgraph, M::Int64, tau::Int64 ; outfile="", T= 100000, beta = 1.0, verbose = false,conf_init=rand(1:graph.q, graph.L), nprocs = 1)

Sample `M` configurations from probability distribution defined by `graph`. Number of MCMC steps between configurations is `tau`. 

Keyword parameters: 
- outfile: Default "". Location to save sample. 
- T: Default 50*tau. Initial number of iterations. 
- beta: Default 1.0. Inverse temperature
- conf_init: Default random. Initial configuration. 
- nprocs: Default 1. Number of parallel MCMC chains.
"""
function doMCMC(graph::DCAgraph, M::Int64, tau::Int64 ; outfile="", T= 50*tau, beta = 1.0, verbose = false,conf_init=rand(1:graph.q, graph.L), nprocs = 1)

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
    # conf_init = rand(1:graph_local.q, graph_local.L)
    conf_end = zeros(Int64, graph_local.L)
    samplefromgraph!(graph_local, conf_init, conf_end, T)
    sample[1,:] = conf_end
    verbose ? println("done!") : print("")

    # Sampling
    verbose ? println("Sampling") : print("")
    @time for m = 1:M-1
        if verbose && mod(m+1,500)==0
            @printf("It %d/%d           \r",m+1,M)
        end
        samplefromgraph!(graph_local, sample[m,:], conf_end, tau)
        sample[m+1,:] .= conf_end
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

Sample for `tau` iterations from probability defined by `g`, starting with configuration `conf_init` and storing final configuration in `conf_end`. 
"""
function samplefromgraph!(g::DCAgraph, conf_init::Array{Int64,1}, conf_end::Array{Int64,1}, tau::Int64)
	E = 0.
	q = g.q
	L = g.L

	copyto!(conf_end, conf_init)
	@fastmath @inbounds for t in 1:tau
		for i in 1:L
			a = conf_end[i]
			b = Random.rand(1:q) 
			while b==a
				b = Random.rand(1:q)
			end
        	id_i_a = (i-1)*q+a
        	id_i_b = (i-1)*q+b
        	E = g.h[id_i_a] - g.h[id_i_b]
        	for j = 1:L
        	    if j!=i
        	        E += g.J[(j-1)*q+conf_end[j], id_i_a] - g.J[(j-1)*q+conf_end[j], id_i_b]
        	    end
        	end

        	if E<=0 || exp(-E) > Random.rand() 
        		conf_end[i] = b
        	end
    	end
	end
end

"""
	estimatetau(g::DCAgraph)

Attempt to estimate reasonable number of iterations between samples. Based on autocorrelation. 
- Conservative: if the autocorrelation of the most autocorrelated spin is smaller than 1/e, eq. is reached. 
- Fast: if the average absolute autocorrelation etc... 
"""
function estimatetau(g::DCAgraph ; itau = 2, M = 1000, threshold = 1/2.7, mode = "conservative")

	t = doMCMC(g, M, itau, T=10000)
	ac = autocorr(t, g.q)
	out = Int64(0)

	if mode == "conservative"
		score = vec(findmax(ac, dims=2)[1])
	elseif mode == "fast"
		score = vec(mean(abs.(ac),dims=2))
	end
	if typeof(findnext(x->x<threshold, score, 1)) != Nothing
		out =  Int64(findnext(x->x<threshold, score, 1)) * itau
	else
		out =  size(ac,1) * itau
	end

	return Int64(out)
end

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

end