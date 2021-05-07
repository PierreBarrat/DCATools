"""
	mutable struct SingleMut

# Summary
Describe a single mutant.

# Fields
- `i::Int64`: Position of the mutation.
- `a::Int64`: New state of the configuration.
- `E::Float64`: Energy of this mutant in a DCA model.
- `fitness::Float64`: Fitness of this mutant.
"""
mutable struct SingleMut
	i::Int64
	a::Int64
	E::Float64
	fitness::Float64
end

"""
	mutable struct Mutant

# Summary
Describe a general mutant.

# Fields
- `smut::Array{SingleMut,1}`: Array of corresponding single mutants.
- `E::Float64`: Energy of the mutant according to a DCA model
- `dE::Float64`: Epistatic effect of energies. Difference between `E` and `+(smut.E)`
- `fitness::Float64`: Fitness of the mutant.
- `dfitness::Float64`: Epistatic effect of fitnesses. Difference between `fitness` and `+(smut.fitness)`
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

Construct empty `Mutant`.
"""
Mutant() = Mutant(Array{SingleMut,1}(undef, 0), 0, 0, 0, 0)

"""
	mutable struct MutData

# Summary
Collection of `Mutant`. Represent a mutagenesis experiment.

# Fields
- `wt::Array{Int64,1}`: Wild-type sequence for this mutagenesis experiment.
- `E_wt::Float64`: Energy of the wt.
- `fitness_wt::Float64`: Fitness of the wt.
- `mutant::Array{Mutant,1}`: Array of `Mutant`.
- `L::Int64`: Length of sequences.
- `q::Int64`: Number of states.
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

Construct empty `MutData`.
"""
MutData() = MutData(Array{Int64,1}(undef, 0),0.,0., Array{Mutant,1}(undef, 0),0,0)
