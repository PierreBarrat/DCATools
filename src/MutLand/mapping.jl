"""
	mapsinglemut!(md::MutData)

For all mutants `X` in `md.mutant`, find corresponding single mutants also in `md.mutant` and store their fitnesses and energies in `X.smut`.
"""
function mapsinglemut!(md::MutData)
	idsmut = findsinglemut(md)[2]
	for (n,mut) in enumerate(md.mutant)
		# If mut is a single mutant, no problem
		if size(mut.smut,1)==1
			mut.smut[1].fitness = mut.fitness
			mut.smut[1].E = mut.E
		else
		# Else, we should look for all single mutants `smutant` in `md` that might correspond to `mut`
			mlist = map(x->(x.i,x.a), mut.smut) # List of single mutations `(i,a)` in `mut`
			for smutant in md.mutant[idsmut]
				if in((smutant.smut[1].i, smutant.smut[1].a), mlist) # smutant.smut[1] --> contains `(i,a)` of this single mutant `smutant`
					idx = findall(in([(smutant.smut[1].i, smutant.smut[1].a)]),mlist)[1] # Find the index in `mut.smut` corresponding to `smut.smut`
					mut.smut[idx].fitness = smutant.fitness
					mut.smut[idx].E = smutant.E # We have data for this single mutant, so energy can be set to 0
					# Note : smutant.smut[1] is not guaranteed to contain a fitness or E value at this point!!
				end
			end
		end
	end
end


# Note
# Energies should be in the wt gauge, ie differences to the wt energy.
# --> I don't care about the energy of the wt in epistasis, I just add differential effects wr/ to the wt
"""
	computeepistasis!(mut::Mutant)

Computes difference between sum of energies and fitnesses of the single mutants composing `mut`, and that of `mut` itself.
"""
function computeepistasis!(mut::Mutant)
	mut.dfitness = mut.fitness - mapreduce(x->x.fitness, +, mut.smut, init=0.)
	mut.dE = mut.E - mapreduce(x->x.E, +, mut.smut, init=0.)
end

"""
	computeepistasis!(md::MutData)

Apply `computeepistasis!` on all mutants in `md`.
"""
function computeepistasis!(md::MutData)
	map(x->computeepistasis!(x), md.mutant)
end

"""
	findsinglemut(md::MutData, i::Int64)

Find all single mutants in `md` corresponding to position `i`.
"""
function findsinglemut(md::MutData, i::Int64)
	slist = Array{Mutant, 1}(undef, 0)
	idlist = Array{Int64,1}(undef, 0)
	for (n,mut) in enumerate(md.mutant)
		if size(mut.smut,1)==1 && mut.smut[1].i == i
			push!(slist, mut)
			push!(idlist, n)
		end
	end
	return (slist,idlist)
end

"""
	findsinglemut(md::MutData)

Find all single mutants in `md`.
"""
function findsinglemut(md::MutData)
	slist = Array{Mutant, 1}(undef, 0)
	idlist = Array{Int64,1}(undef, 0)
	for (n,mut) in enumerate(md.mutant)
		if size(mut.smut,1)==1
			push!(slist, mut)
			push!(idlist, n)
		end
	end
	return (slist,idlist)
end

"""
	finddoublemut(md::MutData, i::Int64, j::Int64)

Find all double mutants in `md` corresponding to position `i` and `j`.
"""
function finddoublemut(md::MutData, i::Int64, j::Int64)
	dlist = Array{Mutant, 1}(undef, 0)
	idlist = Array{Int64, 1}(undef, 0)
	for (n,mut) in enumerate(md.mutant)
		if size(mut.smut,1)==2 && ((mut.smut[1].i, mut.smut[2].i) == (i,j) || (mut.smut[1].i, mut.smut[2].i) == (j,i))
			push!(dlist, mut)
			push!(idlist, n)
		end
	end
	return (dlist,idlist)
end


"""
	finddoublemut(md::MutData)

Find all double mutants in `md`.
"""
function finddoublemut(md::MutData)
	dlist = Array{Mutant, 1}(undef, 0)
	idlist = Array{Int64, 1}(undef, 0)
	for (n,mut) in enumerate(md.mutant)
		if size(mut.smut,1)==2
			push!(dlist, mut)
			push!(idlist, n)
		end
	end
	return (dlist,idlist)
end
