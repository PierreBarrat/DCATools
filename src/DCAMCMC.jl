module DCAMCMC
	
using DCATools
using DelimitedFiles
using Printf
using Statistics

export doMCMC2p, samplefromgraph!, estimatetau

"""
	doMCMC2p(graphfile::String,outfile::String, format::String, q::Int64,M::Int64,t::Int64 ; T= 10000, beta = 1.0, verbose = true, conf_init=rand(1:q, L))

Takes file as input. Untested yet. 
"""
function doMCMC2p(graphfile::String,outfile::String, format::String, q::Int64,M::Int64,t::Int64 ; T= 10000, beta = 1.0, verbose = true, conf_init=rand(1:q, L))
    g = readparam(graphfile, format=format, q=q)
    doMCMC2p(g, M, t, outfile=outfile, T=T, beta=beta, verbose=verbose, conf_init=conf_init)
    return(sample)
end

"""
	doMCMC2p(graph::DCAgraph, M::Int64, tau::Int64 ; outfile="", T= 100000, beta = 1.0, verbose = false,conf_init=rand(1:graph.q, graph.L))

Sample `M` configurations from probability distribution defined by `graph`. Number of MCMC steps between configurations is `tau`. 

Keyword parameters: 
- outfile: Default "". Location to save sample. 
- T: Default 50*tau. Initial number of iterations. 
- beta: Default 1.0. Inverse temperature
- conf_init: Default random. Initial configuration. 
"""
function doMCMC2p(graph::DCAgraph, M::Int64, tau::Int64 ; outfile="", T= 50*tau, beta = 1.0, verbose = false,conf_init=rand(1:graph.q, graph.L))

	# Unnecessary allocation here?  
	graph_local = deepcopy(graph)
    if beta != 1.
        graph_local.J *= beta
        graph_local.h *= beta
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
    for m = 1:M-1
        if verbose && mod(m+1,500)==0
            @printf("It %d/%d           \r",m+1,M)
        end
        samplefromgraph!(graph_local, sample[m,:], conf_end, tau)
        sample[m+1,:] = conf_end
    end
    verbose ? println("") : print("")
    verbose ? println("Done") : print("")

    # Writing output
    if outfile!=""
        if outformat=="simple"
            writedlm(outfile,sample,' ')
        elseif outformat=="Matteo"
            ff = open(outfile,"w")
            write(ff,@sprintf("%d %d %d\n",M,graph_local.L,graph_local.q))
            writedlm(ff, sample-1, ' ')
            close(ff)
        end
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
		i = rand(1:L)
		a = conf_end[i]
		b = rand(1:q)
		while b==a
			b = rand(1:q)
		end
        id_i_a = (i-1)*q+a
        id_i_b = (i-1)*q+b
        E = g.h[id_i_a] - g.h[id_i_b]
        for j = 1:L
            if j!=i
                E += g.J[(j-1)*q+conf_end[j], id_i_a] - g.J[(j-1)*q+conf_end[j], id_i_b]
            end
        end

        if E<=0 || exp(-E) > rand(Float64) 
        	conf_end[i] = b
        end
	end
end

"""
	estimatetau(g::DCAgraph)

Attempt to estimate reasonable number of iterations between samples. Based on autocorrelation. 
- Conservative: if the autocorrelation of the most autocorrelated spin is smaller than 0.1, eq. is reached. 
- Fast: if the average absolute autocorrelation etc... 
"""
function estimatetau(g::DCAgraph ; itau = 50, M = 5000, threshold = 0.1, mode = "conservative")
	t = doMCMC2p(g, M, itau, T=10000)
	uptau = itau

	ac = autocorr(t, g.q)

	if mode == "conservative"
		score = vec(findmax(ac, dims=2)[1])
	elseif mode == "fast"
		score = vec(mean(abs.(ac),dims=2))
	end
	if typeof(findnext(x->x<threshold, score, 1)) != Nothing
		return findnext(x->x<threshold, score, 1) * itau
	else
		return size(ac,1) * itau
	end
end

"""
"""
function hdist(s1::Array{Int64,1}, s2::Array{Int64,1})
	out = 0.
	for i in 1:size(s1,1)
		out += Int64(s1[i]!=s2[i])
	end
	return out/size(s1,1)
end

"""
"""
function autocorr(sample::Array{Int64,2}, q::Int64)
	(M,L) = size(sample)
	navmin = Int64(M-round(M/10))
	f1 = computefreqs(sample)[1]

	ac = zeros(Float64, M-navmin, L)
	for t in 1:(M-navmin)
		for m in 1:(M-t)
			for i in 1:L
				ac[t,i] += Float64(sample[m,i] == sample[m+t,i]) 
			end
		end
		ac[t,:] /= (M-t)
		ac[t,:] -= vec(sum(reshape(f1,21,31).^2, dims=1)')
	end
	return ac
end

end