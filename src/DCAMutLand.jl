module DCAMutLand


export SingleMut, Mutant, MutData, readmutdata, mapsinglemut!, findsinglemut, finddoublemut, computeepistasis!

using DCATools

"""
	SingleMut

# Summary 
Describe a single mutant. 
# Fields
- i::Int64 - Position of the mutation.
- a::Int64 - New state of the configuration. 
- E::Float64 - Energy of this mutant in a DCA model.
- fitness::Float64 - Fitness of this mutant. 
"""
mutable struct SingleMut
	i::Int64
	a::Int64
	E::Float64
	fitness::Float64
end

"""
	Mutant
# Summary
Describe a general mutant.
# Fields
- smut::Array{SingleMut,1} - Array of corresponding single mutants. 
- E::Float64 - Energy of the mutant according to a DCA model
- dE::Float64 - Epistatic effect of energies. Difference between `E` and `+(smut.E)`
- fitness::Float64 - Fitness of the mutant.
- dfitness::Float64 - Epistatic effect of fitnesses. Difference between `fitness` and `+(smut.fitness)`
"""
mutable struct Mutant
	smut::Array{SingleMut,1}
	E::Float64
	dE::Float64 
	fitness::Float64
	dfitness::Float64
end
"""
	Mutant()

Basic constructor for `Mutant`. Constructs empty `Mutant`.
"""
Mutant() = Mutant(Array{SingleMut,1}(undef, 0), 0, 0, 0, 0)

"""
	MutData

# Summary
Collection of `Mutant`. Represent a mutagenesis experiment.
# Fields
- wt::Array{Int64,1} - Wild-type sequence for this mutagenesis experiment. 
- E_wt::Float64 - Energy of the wt. 
- fitness_wt::Float64 - Fitness of the wt.
- mutant::Array{Mutant,1} - Array of `Mutant`. 
- L::Int64 - Length of sequences.
- q::Int64 - Number of states. 
"""
mutable struct MutData
	wt::Array{Int64,1}
	E_wt::Float64
	fitness_wt::Float64
	mutant::Array{Mutant,1}
	L::Int64
	q::Int64
end

"""
	MutData()

Basic constructor for `MutData`. 
"""
MutData() = MutData(Array{Int64,1}(undef, 0),0.,0., Array{Mutant,1}(0),0,0)


"""
	readmutdata(infile::String ; mapsingle = false)

Read data stored in `infile` to a `MutData` structure. Mapping single mutants `SingleMutant` of each `Mutant` to their value found in the data is optional. 

Note that energies and epistatic effect are not computed here. 
"""
function readmutdata(infile::String ; mapsingle = false)

	lines = open(infile) do f 
		lines = readlines(f)
	end

	mutdata = MutData()

	for (cl,l) in enumerate(lines)
		# Comments #
		if l=="" || l[1]=='#' 

		# Info (wt, L, q) #
		elseif l[1]=='>'
			if l[2:3]=="wt"
				mutdata.wt = Int64.(parse.(split(l," ")[2:end]))
				# println("wild type: $(mutdata.wt)")
			elseif l[2:5]=="info"
				mutdata.L = Int64(parse(split(l," ")[2]))
				mutdata.q = Int64(parse(split(l," ")[3]))
				# println("L = $(mutdata.L), q = $(mutdata.q)")
			else
				println("MutationalData.jl - readmutdata: line $cl begins with '>' but does not contain any known identifier.")
			end

		# Data # 
		else 
			cmut = Mutant()
			pl = parse.(split(l, " "))
			N = Int64(pl[1])
			for n in 1:N
				push!(cmut.smut, SingleMut(Int64(pl[2*n]), Int64(pl[2*n+1]), NaN, NaN)) 
				# Uninitialized single mutants are set to `NaN`. If they are never observed, this will propagate to the `epistasis` value of `cmut`
				# --> using NaN as missing value (for Float)
			end
			cmut.fitness = pl[end]
			push!(mutdata.mutant, cmut)
		end
	end

	# At this stage, all mutants are stored, but fitnesses of single mutations for multiple mutants are not determined
	# If `mapsingle` is true, mapping between single and multiple mutants will be made. If energies are not yet computed, one might have to do this twice. 
	if mapsingle
		mapsinglemut!(mutdata)
	end
	return mutdata
end


"""
	function mapsinglemut!(md::MutData)

For all mutants `X` in `MutData.mutant`, find corresponding single mutants also in `MutData.mutant` and store their fitnesses and energies in `X.smut`. 
"""
function mapsinglemut!(md::MutData)
	idsmut = findsinglemut(md)[2]
	for (n,mut) in enumerate(md.mutant)
		# If mut is a single mutant, no problem
		if size(mut.smut,1)==1
			mut.smut[1].fitness = mut.fitness
			mut.smut[1].E = mut.E
		# Else, we should look for all single mutants `smut` in `md` that might correspond to `mut`
		else
			mlist = map(x->(x.i,x.a), mut.smut) # List of single mutations `(i,a)` in `mut` 
			for smut in md.mutant[idsmut]
				if in((smut.smut[1].i, smut.smut[1].a), mlist) # smut.smut[1] --> contains `(i,a)` of this single mutant `smut`
					idx = findin(mlist,[(smut.smut[1].i, smut.smut[1].a)])[1] # Find the index in `mut.smut` corresponding to `smut.smut`
					mut.smut[idx].fitness = smut.fitness 
					mut.smut[idx].E = smut.E # We have data for this single mutant, so energy can be set to 0 
					# Note : smut.smut[1] is not guaranteed to contain a fitness or E value at this point!! --> use smut.XXX
				end
			end
		end
	end
end


# Note
# Energies should be in the wt gauge, ie differences to the wt energy. 
# --> I don't care about the energy of the wt in epistasis, I just add differential effects wr/ to the wt
"""
	function computeepistasis!(mut::Mutant)

Computes difference between sum of energies and fitnesses of the single mutants composing `mut`, and that of `mut` itself.
"""
function computeepistasis!(mut::Mutant)
	mut.dfitness = mut.fitness - mapreduce(x->x.fitness, +, 0., mut.smut)
	mut.dE = mut.E - mapreduce(x->x.E, +, 0., mut.smut)
end


"""
	function findsinglemut(md::MutData, i::Int64)

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
	function findsinglemut(md::MutData)

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
	function finddoublemut(md::MutData, i::Int64, j::Int64)

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
	function finddoublemut(md::MutData)

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


include("DCAMutLand/modeltools.jl")

end