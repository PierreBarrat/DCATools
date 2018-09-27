## Basics
DCA parameters could be stored in a structure. 
```
struct graph
    J::Array{Float64,2}
    h::Array{Float64,1}
    L::Int64
    q::Int64
end
```
This makes it a standalone object, with all the information needed to manipulate it easily. 
*We have to agree on a representation of couplings/fields : 4D or 2D?*

## Input / Output
#### Reading DCA parameters from a file -- ok
General form of the function
`function readparam(infile::String ; format = XXX)`
Different formats for the input file should be handled: 

- "Matteo's" format, style `J i j a b val` and `h i a val`
- A matrix format, with a line of the file representing the vector `J[i][:][a][:]` (in a 4D notation to make it clear). The last line represents the fields. 
- others? 
- 
*Pros and cons of each format?* 

The return of `readparam` should be of type `Graph`. 
Depending on the format, it could be necessary to indicate the value of `q` or `L` to the function. A fields vector of size 10 could be `L=5, q=2` or `L=2, q=5` ... 

#### Writing DCA parameters to a file -- ok
`function writeparam(outfile::String, g::Graph ; format = XXX)`
Same formats as above. 

#### Reading alignments -- half
Maybe different functions, depending whether the input is already in a numerical format, or still in a fasta format with letters. 

- `function readmsanum(file::String ; format=XXX, header=true/false)` -- ok

with format:
- starting at zero
- starting at one
Matteo's BM uses a header to indicate number/length of sequences, and number of states
*Need format where sequences are numerical, but the `>` comment line of fasta is still here?* 

- `function readfasta(file::String ; mapping=XXX)`
I guess we would all be better off using the **same** mapping from amino acids to numbers.
*By the way, we don't all deal with proteins. How should we handle this?*
*We could add an alphabet string as optional input*


## Utilities
#### Pairwise frequencies for alignment - OK
This should have several formats, for different inputs. 

- Alignment is in a file, in numerical format.
- Alignment is in a file, in fasta format.
- Alignment is in an `Array{Int64,2}` already loaded into julia.

And additionally, options for weights.

- Weights are in a file of `M` lines.
- Weights are in an `Array{Float64,1}` already loaded into julia.
- No weights should be used.
- Weights have to be computed.

And finally in this last scenario:

- Weights should be computed and stored in a file.
- Weights should be computed and then discarded.

Putting aside the case where the msa is in a fasta file: 

**Base routine, where computation is carried:**

`function computefreqs(Y::Array{Int64,2}, w::Array{Float64,1}, q::Int64)`

**Then, different functions to call this one.**

1. `function computefreqs(Y::Array{Int64,2} ; q = 21, reweighting=false, weights=[])`. 
Handles different cases
 - If `reweighting=true`, calls a routine to compute weights. If `weights` is a `String`, tries to save weights there.
 - If `reweighting=false` and `weights` is a `String`, tries to read weights there. 
 - If `reweighting=false` and `weights` is an `Array{Float64,1}`, calls the base routine directly. 
2. `function computefreqs(infile::String ; q = 21, reweighting=false, weights=[])`

#### Switching gauge
One keyword argument for each gauge. Let's start with 0-sum gauge since it's easy and commonly used for contact prediction. 
`function switchgauge(g::DCAgraph ; gauge="0sum")`
Return a `DCAgraph`. 
In place change with `switchgauge!` could easily be implemented too. 

#### Computing energies
Given an alignment `Array{Int64,2}` or a single sequence `Array{Int64,1}`, compute corresponding energies with a given dca model. 
`function computeenergies(g::DCAgraph, Y::Array{Int64,2})`
`function computeenergies(g::DCAgraph, Y::Array{Int64,1})`
Return an `Array{Float64,1}` containing said energies.
