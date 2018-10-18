module DCAMCMC
	
using DCATools
using DelimitedFiles
using Printf

export doMCMC2p, samplefromgraph! 

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




end