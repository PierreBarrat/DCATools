#=
# IMPORTANT NOTE
# All energies should be computed with respect to the wild-type (*ie* wt should have null energy)
# If not, computing epistasis for an arbitrary mutant becomes a mess
# However, this means that data for different wt cannot be compared: they are not in the same gauge!!
# This is in line with how deltaG are computed in e.g. rama's paper: wr/ to the wt.
# --> if deltaG are averaged over homologs, energy changes wr/ to wt can be averaged as well! However, they are not indicative of a probability since they do not come from the same gauge
#
# The other option would have been to store energies in a 0-sum gauge, and also the energy change wr/ to the wt.
# This creates additional fields in both singlemut and mutant ...
=#


"""
	computeenergies!(md::MutData, g::DCAGraph)

Compute energies of all `Mutant` in `md.mutant` using `g`.
"""
function computeenergies!(md::MutData, g::DCAGraph)
	md.E_wt = computeenergies(g, md.wt)
	for mut in md.mutant
		computeenergies!(mut, g, md.wt, md.E_wt)
	end

	return nothing
end

"""
	computeenergies!(mut::Mutant, g::DCAGraph, wt::Array{Int64,1}, E_wt::Float64)

Compute energies of `mut` and of all single mutants in `mut.smut`. Compute epistatic effects as well using `E_wt`.
"""
function computeenergies!(mut::Mutant, g::DCAGraph, wt::Array{Int64,1}, E_wt::Float64)
	mseq = copy(wt)
	mseqt = copy(wt)
	for smut in mut.smut
		mseq[smut.i] = smut.a
		mseqt[smut.i] = smut.a
		smut.E = energy(g, mseq) .- E_wt
		mseq[smut.i] = wt[smut.i]
	end
	mut.E = energy(g, mseqt) .- E_wt

	return nothing
end

"""
	mapenergies!(md::MutData, g::DCAGraph)

Map energies values of mutants to fitness values. Details of the mapping can be found in
	Coevolutionary landscape inference and the context-dependence of mutations in beta-lactamase TEM-1
	Matteo Figliuzzi, Hervé Jacquier, Schug Alexander, Olivier Tenaillon, Martin Weigt
Energies in `md` are modified in the process.

Output `mapping` is a dictionary such that `mapping[E] = fitness`.
"""
function mapenergies!(md::MutData, g::DCAGraph)
	computeenergies!(md, g)
	mapping = Dict{Float64, Float64}()
	fitlist = Array{Float64,1}(undef, 0)
	Elist = Array{Float64,1}(undef, 0)
	for m in md.mutant
		push!(fitlist, m.fitness)
		push!(Elist, m.E)
	end
	sort!(fitlist)
	sort!(Elist,rev=true)
	map(x->mapping[Elist[x]]=fitlist[x], 1:size(fitlist,1))
	for m in md.mutant
		m.E = mapping[m.E]
	end

	return mapping
end

"""
	mapenergies(md::MutData, g::DCAGraph)

Map energies values of mutants to fitness values. Details of the mapping can be found in
	Coevolutionary landscape inference and the context-dependence of mutations in beta-lactamase TEM-1
	Matteo Figliuzzi, Hervé Jacquier, Schug Alexander, Olivier Tenaillon, Martin Weigt

Output `mapping` is a dictionary such that `mapping[E] = fitness`.
"""
function mapenergies(md::MutData, g::DCAGraph)
	md_ = deepcopy(md)
	computeenergies!(md_, g)
	mapping = Dict{Float64, Float64}()
	fitlist = Array{Float64,1}(undef, 0)
	Elist = Array{Float64,1}(undef, 0)
	for m in md_.mutant
		push!(fitlist, m.fitness)
		push!(Elist, m.E)
	end
	sort!(fitlist)
	sort!(Elist,rev=true)
	map(x->mapping[Elist[x]]=fitlist[x], 1:size(fitlist,1))

	return mapping
end
