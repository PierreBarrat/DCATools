export bmstep!, bmlearn

let verbose::Bool = false
	global v() = verbose
	global vv() = vverbose
	global set_verbose(v) = (verbose = v)
	global set_vverbose(v) = (vverbose = v)
end

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
"""
function writeinfo(infofile::String, meta::BMmeta, ginit::Bool, mutants::Bool)
	open(infofile, "w") do of
		write(of, "Parameters of the bmlearn run:\n")
		for f in fieldnames(BMmeta)
			if f != :comment
				write(of, "$f: $(getproperty(meta, f))\n")
			end
		end
		write(of, "Initial non-zero graph provided: $ginit\n")
		write(of, "Mutational data provided: $mutants\n")
		write(of, "\nComment: $(meta.comment)")
	end
	nothing
end

"""
	bmlearn(f1::Array{Float64,1}, f2::Array{Float64,2}, L::Int64, q::Int64;
			ginit::DCAgraph = DCAgraph(L,q), gradinit=DCAgrad(L, q), mutants::MutData = MutData(),
			kwargs...)

Learn a Boltzmann machine with target frequencies `f1` and `f2`. `L` and `q` are the length of sequences and number of characters of the alphabet. 
`ginit` provides a starting point for the BM. If provided, `mutants::MutData` will be used to integrate mutational information (*e.g.* from a deep mutational scan). 

Parameters guiding the learning process are provided as additional keyword arguments and passed to `BMmeta`. See `?BMmeta` for help. 
"""
function bmlearn(f1::Array{Float64,1}, f2::Array{Float64,2}, L::Int64, q::Int64;
			ginit::DCAgraph = DCAgraph(L,q), gradinit=DCAgrad(L, q), mutants::MutData = MutData(),
			kwargs...)
	try 
		BMmeta(;kwargs...)
	catch err
		@warn "Unsupported keyword argument. See `?BMmeta` for help."
		error(err)
	end
	bmlearn(f1, f2, L, q, BMmeta(;kwargs...), ginit=ginit, gradinit=gradinit, mutants=mutants)
end

"""
	bmlearn(f1::Array{Float64,1}, f2::Array{Float64,2}, L::Int64, q::Int64, meta::BMmeta = BMmeta();
			ginit::DCAgraph = DCAgraph(L,q), gradinit=DCAgrad(L, q), mutants::MutData = MutData())
"""
function bmlearn(f1::Array{Float64,1}, f2::Array{Float64,2}, L::Int64, q::Int64, meta::BMmeta;
	ginit::DCAgraph = DCAgraph(L,q), gradinit=DCAgrad(L, q), mutants::MutData = MutData())
	
	# Initializing save directory
	mkpath(meta.savefolder)
	logfile = "$(meta.savefolder)/$(meta.logfile)"

	# Writing info about run
	writeinfo("$(meta.savefolder)/$(meta.infofile)", meta, ginit!=DCAgraph(L,q), !isempty(mutants.mutant))

	# Setting verbosity
	set_verbose(meta.verbose)


	# Initializing meta data and log
	samplesize = meta.Minit
	bmlog = BMlog()
	bmlog.samplesize = samplesize

	# Initializing gradient and sample
	g = deepcopy(ginit)	
	cgrad = deepcopy(gradinit)
	sample = bminit!(cgrad, g, f1, f2, meta, bmlog)

	# Copying mutants
	md = deepcopy(mutants)
	# if !isempty(md.mutant)
	# 	switchgauge!(g, gauge="wt", wt = md.wt)
	# end

	# Header for meta.logfile
	writelog(meta.logfile)

	for it in 1:meta.nit
		if mod(it,meta.update_tau) == 1
			bmlog.tau = 0
		end
		cgrad = bmstep!(g, f1, f2, md, cgrad, meta, bmlog)
		v() && println(" --- It. $it out of $(meta.nit) ---")
		v() && println("Norm of gradient: $(bmlog.gradnorm)")
		writelog(logfile, bmlog)
		if mod(it, meta.saveparam)==0 && meta.savefolder!=""
			writeparam("$(meta.savefolder)/DCABM_it$(it)_mat.txt", g, format="mat")
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
		tau = min(estimatetau(g)[1], g.L)
	end
	return doMCMC(g, M, tau, T = 10*tau, nprocs = nprocs), tau
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
	if !ismissing(meta.l2)
		reg = computel2(g, meta.l2)
	else
		reg = computel2(g, meta.l2J, meta.l2h)
	end
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

