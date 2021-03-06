export bmstep!, bmlearn

"""
"""
function writelog(logname::String, bmlog::BMlog)
	ff = open(logname,"a")
	outstring = @sprintf("%d    %d    %.4f    %.4f    %.4f    %.4f    %.4f    %.4f    %.4f    %.4f    %.4f\n",bmlog.samplesize,bmlog.tau, bmlog.gradnorm, bmlog.gradnormh, bmlog.gradnormJ, bmlog.gradconsth, bmlog.gradconstJ, bmlog.corcor, bmlog.slopecor, bmlog.cormag, bmlog.cormutants)
	write(ff, outstring)
	# write(ff, "$(bmlog.samplesize)\t$(bmlog.tau)\t$(bmlog.gradnorm)\t$(bmlog.gradnormh)\t$(bmlog.gradnormJ)\t$(bmlog.gradconsth)\t$(bmlog.gradconstJ)")
	# write(ff, "\t$(bmlog.corcor)\t$(bmlog.slopecor)\t$(bmlog.cormag)")
	# write(ff, "\n")
	close(ff)
end
"""
"""
function writelog(logname::String)
	ff = open(logname,"w")
	write(ff, "M    tau    gradient_norm    gradnormh    gradnormJ    gradient_consistency_h    gradient_consistency_J")
	write(ff, "    corcor    slopecor    cormag    cormutants")
	write(ff, "\n")
	close(ff)
end


"""
	bmlearn(f1::Array{Float64,1}, f2::Array{Float64,2}, L::Int64, q::Int64)

Main BM learning function. Keywords correspond to starting point and meta parameters of the learning process.

Some useful ones
- `nit`: Number of iterations of the gradient descent. 
- `savefolder`: folder for output files (log and parameters). If unspecified or `""`, current directory is used.
- `saveparam`: number of iterations after which parameters are saved.  
- `logfile`: name of logfile. 
- `ginit`: `DCAgraph` object. Initial values of parameters. 
- `gradinit`: `DCAgrad` object. Initial value of gradient. Useful to restart a learning. Can also be used to train the model only on a subset of parameters by setting initial `stepJ` and `steph` to non zero values for those only. 
"""
function bmlearn(f1::Array{Float64,1}, f2::Array{Float64,2}, L::Int64, q::Int64 ; 
	ginit::DCAgraph = DCAgraph(L,q), gradinit=DCAgrad(L, q), 
	l2 = 0.01, l1 = 0.,
	samplesize = 1000, 
	basestepJ = 0.05, basesteph = 0.05,  stepJmax = 0.5, stephmax = 0.5,
	aJup = 1.2, aJdown = 0.5, ahup = 1.2, ahdown = 0.5, 
	adaptMup = 1.05, Mmax = 100000, 
	update_tau = 10,
	saveparam = 10, savefolder="",
	mutants::MutData = MutData(), integrative_lambda=0, integrative_M = 1, 
	nit = 50, 
	nprocs = 1,
	logfile="log.txt")
	
	# Initializing save directory
	if savefolder!=""
		mkpath(savefolder)
		logfile = "$(savefolder)/$(logfile)"
	end

	# Initializing meta data and log
	meta = BMmeta(l2, l1, basestepJ, basesteph, stepJmax, stephmax, aJup, aJdown, ahup, ahdown, adaptMup, Mmax, integrative_lambda, integrative_M, saveparam, nprocs)
	bmlog = BMlog()
	bmlog.samplesize = samplesize

	# Initializing gradient and sample
	g = deepcopy(ginit)	
	cgrad = deepcopy(gradinit)
	sample = bminit!(cgrad, g, f1, f2, meta, bmlog)

	# Copying mutants
	md = deepcopy(mutants)
	if !isempty(md.mutant)
		switchgauge!(g, gauge="wt", wt = md.wt)
	end

	# Header for logfile
	writelog(logfile)

	for it in 1:nit
		if mod(it,update_tau) == 1
			bmlog.tau = 0
		end
		cgrad = bmstep!(g, f1, f2, md, cgrad, meta, bmlog)
		println(" --- It. $it out of $nit ---")
		println("Norm of gradient: $(bmlog.gradnorm)")
		writelog(logfile, bmlog)
		if mod(it, meta.saveparam)==0 && savefolder!=""
			writeparam("$(savefolder)/DCABM_it$(it)_mat.txt", g, format="mat")
		elseif mod(it, meta.saveparam)==0
			writeparam("DCABM_it$(it)_mat.txt", g, format="mat")
		end
	end
	return g, cgrad, bmlog
end



"""
	bmsample(g::DCAgraph, M::Int64)
"""
function bmsample(g::DCAgraph, M::Int64, nprocs, tau::Int64)
	if tau == 0
		tau = 3*estimatetau(g, mode="fast", nprocs=nprocs)
	end
	return doMCMC(g, M, tau, T = 20*tau, nprocs = nprocs), tau
end

"""
	bmstep!(g::DCAgraph, f1::Array{Float64,1}, f2::Array{Float64,2}, md::MutData, prevgrad::DCAgrad, meta::BMmeta, bmlog::BMlog)

Compute gradient for current graph `g` and updates it. Target frequencies are `f1` and `f2`. Local mutational data is `md`. Return computed gradient. 
"""
function bmstep!(g::DCAgraph, f1::Array{Float64,1}, f2::Array{Float64,2}, md::MutData, prevgrad::DCAgrad, meta::BMmeta, bmlog::BMlog)

	# Compute new sample
	sample, tau = bmsample(g, bmlog.samplesize, meta.nprocs, bmlog.tau)
	bmlog.tau = tau


	# Compute gradient from frequency difference and l2 regularization
	freqgrad, p1, p2 = computegradient(sample, f1, f2, g.q)
	reg = computel2(g, meta.l2)
	gradtot = freqgrad + reg

	# If l1 regularization exists, add it to gradient
	if meta.l1 != 0
		gradtot += computel1!(g, gradtot, meta.l1)
	end

	# If we're also fitting mutants, add it to gradient
	bmlog.cormutants = 0.
	if !isempty(md.mutant)
		computeenergies!(md, g)
		mapping = mapenergies(md, g)
		bmlog.cormutants = cor(map(x->x.fitness, md.mutant), map(x->mapping[x.E], md.mutant))
		mutgrad = computegradient(md, mapping, meta)
		gradtot += mutgrad
		# We want to stay in the wt gauge --> consider only gradient on non-wt positions
		for i in 1:g.L
			gradtot.gradh[(i-1)*g.q + md.wt[i]] = 0.
			for j in (i+1):g.L
				gradtot.gradJ[(i-1)*g.q + md.wt[i], (j-1)*g.q + md.wt[j]] = 0.
				gradtot.gradJ[(j-1)*g.q + md.wt[j], (i-1)*g.q + md.wt[i]] = 0.
			end
		end
	end

	# Determine step size -- adaptive
	computestepsize!(gradtot, prevgrad, meta)

	# Update parameters
	updateparameters!(g, gradtot)

	# Updating M and computing gradient norm 
	bmlog.gradnorm, bmlog.gradnormh, bmlog.gradnormJ = gradnorm(gradtot)
	bmlog.gradconsth, bmlog.gradconstJ = gradconst(prevgrad, gradtot)
	# updateM!(bmlog, prevgrad, gradtot, meta)
	updateM!(bmlog, meta)

	# Fitting quality at this iteration
	fq = fitquality(f2, f1, p2, p1, g.q)
	bmlog.corcor = fq[1]
	bmlog.slopecor = fq[2]
	bmlog.cormag = fq[4]


	return gradtot
end

"""
	bminit!(grad::DCAgrad, g::DCAgraph, f1::Array{Float64,1}, f2::Array{Float64,2}, meta::BMmeta)

Initialize gradient and sample with null values. Step sizes for gradient are initialized using `meta`. If the gradient was a non-null one, *ie* if gradient values or stepsizes are different from 0, then the full gradient object remains untouched. 
"""
function bminit!(grad::DCAgrad, g::DCAgraph, f1::Array{Float64,1}, f2::Array{Float64,2}, meta::BMmeta, bmlog::BMlog)
	sampleinit = zeros(Int64, bmlog.samplesize, g.L)

	if gradequal(DCAgrad(g.L, g.q), grad)
		grad.stepJ .= meta.basestepJ
		grad.steph .= meta.basesteph
	end
	for i in 1:g.L
		grad.stepJ[(i-1)*g.q .+ (1:g.q), (i-1)*g.q .+ (1:g.q)] .= 0
	end
	return sampleinit
end

"""
	updateM!(bmlog::BMlog, prevgrad::DCAgrad, newgrad::DCAgrad, meta::BMmeta)

Update `M` based on the norM of the gradient. If gradient norm (for `J` only) increased in the last iteration, `M` is increased.
"""
function updateM!(bmlog::BMlog, prevgrad::DCAgrad, newgrad::DCAgrad, meta::BMmeta)
	newgrad_norm = sum(newgrad.gradJ.^2)
	prevgrad_norm = sum(prevgrad.gradJ.^2)
	if newgrad_norm > prevgrad_norm && prevgrad_norm!=0.
		bmlog.samplesize = min(meta.Mmax, Int64(round(meta.adaptMup * bmlog.samplesize)))
	end
end

"""
	updateM!(bmlog::BMlog, meta::BMmeta ; threshold = 0.4)

Update `M` based on gradient consistency. If cosine between the two gradient is smaller than `threshold`, `M` is increased. 
"""
function updateM!(bmlog::BMlog, meta::BMmeta ; threshold = 0.4)
	if bmlog.gradconstJ < threshold
		bmlog.samplesize = min(meta.Mmax, Int64(round(meta.adaptMup * bmlog.samplesize)))
	end
end

