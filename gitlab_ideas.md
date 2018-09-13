### Basics
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

### Input / Output
#### Reading DCA parameters from a file
General form of the function
`function readparam(infile::String ; format = XXX)`
Different formats for the input file should be handled: 
- "Matteo's" format, style `J i j a b val` and `h i a val`
- A matrix format, with a line of the file representing the vector `J[i][:][a][:]` (in a 4D notation to make it clear). The last line represents the fields. 
- others? 
*Pros and cons of each format?* 

The return of `readparam` should be of type `Graph`. 
Depending on the format, it could be necessary to indicate the value of `q` or `L` to the function. A fields vector of size 10 could be `L=5, q=2` or `L=2, q=5` ... 

#### Writing DCA parameters to a file
`function writeparam(outfile::String, g::Graph ; format = XXX)`
Same formats as above. 

#### Reading alignments
Maybe different functions, depending whether the input is already in a numerical format, or still in a fasta format with letters. 

`function readmsanum(file::String ; format=XXX, header=true/false)`
with format:
- starting at zero
- starting at one
Matteo's BM uses a header to indicate number/length of sequences, and number of states
*Need format where sequences are numerical, but the `>` comment line of fasta is still here?* 

`function readfasta(file::String ; mapping=XXX)`
I guess we would all be better off using the **same** mapping from amino acids to numbers.
*By the way, we don't all deal with proteins. How should we handle this?*

