"""
	readmutdata(infile::AbstractString ; mapsingle = false)

Read data stored in `infile` to a `MutData` structure. Mapping single mutants `SingleMutant` of each `Mutant` to their value found in the data is optional.

`infile` should look like this:
```
>wt x1 x2 x3 ... *numerical wt sequence, space separated*
>info L q
n i a j b ... val *where n is the degree of the mutant*
.              .
.              .
.              .
```

Does **not** compute energies or epistatic effects.
"""
function readmutdata(infile::AbstractString ; mapsingle = false)

	lines = open(infile) do f
		lines = readlines(f)
	end

	mutdata = MutData()

	for (cl,l) in enumerate(lines)
		if l[1]=='>' # Info (wt, L, q)
			if l[2:3]=="wt"
				mutdata.wt = Int64.(Meta.parse.(split(l," ")[2:end]))
			elseif l[2:5]=="info"
				mutdata.L = Int64(Meta.parse(split(l," ")[2]))
				mutdata.q = Int64(Meta.parse(split(l," ")[3]))
			else
				println("MutationalData.jl - readmutdata: line $cl begins with '>' but does not contain any known identifier.")
			end
		elseif l!="" && l[1]!='#' # Data
			cmut = Mutant()
			pl = Meta.parse.(split(l))
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


# function export_as_dataframe(md::MutData)

# end
