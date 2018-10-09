
<a id='DCATools.jl-Documentation-1'></a>

# DCATools.jl Documentation

- [DCATools.jl Documentation](index.md#DCATools.jl-Documentation-1)
    - [Data structure](index.md#Data-structure-1)
    - [Alignment tools – alignmenttools.jl](index.md#Alignment-tools-–-alignmenttools.jl-1)
        - [Compute pairwise frequencies.](index.md#Compute-pairwise-frequencies.-1)
        - [Compute weights of sequences.](index.md#Compute-weights-of-sequences.-1)
    - [Input / output - inputoutput.jl](index.md#Input-/-output-inputoutput.jl-1)
        - [Read DCA parameters from file.](index.md#Read-DCA-parameters-from-file.-1)
        - [Write DCA parameters to file.](index.md#Write-DCA-parameters-to-file.-1)
        - [Read a numerical alignment from a file.](index.md#Read-a-numerical-alignment-from-a-file.-1)
    - [Model tools – modeltools.jl](index.md#Model-tools-–-modeltools.jl-1)
        - [Switch gauge of DCA parameters.](index.md#Switch-gauge-of-DCA-parameters.-1)
        - [Compute energies of sequences for given DCA parameters.](index.md#Compute-energies-of-sequences-for-given-DCA-parameters.-1)
        - [Infer a profile model from frequencies or an alignment.](index.md#Infer-a-profile-model-from-frequencies-or-an-alignment.-1)
        - [Compute pseudo-likelihood of data](index.md#Compute-pseudo-likelihood-of-data-1)
    - [Contact prediction – contactprediction.jl](index.md#Contact-prediction-–-contactprediction.jl-1)
        - [Compute Fapc score for given coupling matrix.](index.md#Compute-Fapc-score-for-given-coupling-matrix.-1)
        - [Compute PPV for given scores and distances.](index.md#Compute-PPV-for-given-scores-and-distances.-1)
    - [Miscellaneous – misc.jl](index.md#Miscellaneous-–-misc.jl-1)
        - [Fitting quality of DCA model](index.md#Fitting-quality-of-DCA-model-1)
        - [Three points correlation](index.md#Three-points-correlation-1)


<a id='Data-structure-1'></a>

## Data structure

<a id='DCATools.DCAgraph' href='#DCATools.DCAgraph'>#</a>
**`DCATools.DCAgraph`** &mdash; *Type*.



```
struct DCAgraph
```


<a target='_blank' href='https://github.com/PierreBarrat/DCATools/blob/7b76e84e35cb66ec989505c4dced460cd09285ea/src/DCATools.jl#L7-L9' class='documenter-source'>source</a><br>

<a id='Base.:*' href='#Base.:*'>#</a>
**`Base.:*`** &mdash; *Function*.



```
function *(B, g::DCAgraph)
```

Multiply fields and couplings in `g` by scalar `B`. Useful to change temperature. 


<a target='_blank' href='https://github.com/PierreBarrat/DCATools/blob/7b76e84e35cb66ec989505c4dced460cd09285ea/src/DCATools.jl#L18-L22' class='documenter-source'>source</a><br>


```
function *(B, g::DCAgraph)
```

Multiply fields and couplings in `g` by scalar `B`. Useful to change temperature. 


<a target='_blank' href='https://github.com/PierreBarrat/DCATools/blob/7b76e84e35cb66ec989505c4dced460cd09285ea/src/DCATools.jl#L27-L31' class='documenter-source'>source</a><br>


<a id='Alignment-tools-–-alignmenttools.jl-1'></a>

## Alignment tools – alignmenttools.jl


<a id='Compute-pairwise-frequencies.-1'></a>

### Compute pairwise frequencies.

<a id='DCATools.computefreqs' href='#DCATools.computefreqs'>#</a>
**`DCATools.computefreqs`** &mdash; *Function*.



```
computefreqs(Y::Array{Int64,2}, w::Array{Float64,1}, q::Int64)
```

Base routine for computing frequencies in sample `Y`. `w` is an array containing the weights. 


<a target='_blank' href='https://github.com/PierreBarrat/DCATools/blob/7b76e84e35cb66ec989505c4dced460cd09285ea/src/alignmenttools.jl#L3-L7' class='documenter-source'>source</a><br>


```
computefreqs(Y::Array{Int64,2}; q = findmax(Y)[1], computew=false, weights=[], theta=0.2, saveweights="")
```

Compute pairwise frequencies for an array input. `Y` is an array of `Int64`. Return frequencies `f1` and `f2` and weights.

Keywords: 

  * `q`: default `findmax(Y)[1]`
  * `weights`: default `[]`. If it is a `String`, phylogenetic weights are read from the corresponding file. If it is an `Array{Float64,1}`, they are used directly.
  * `computew`: default `false`. If true, phylogenetic weights are computed, calling the appropriate `computeweights`. `weights` is then ignored.
  * `saveweights` and `theta`: see `computeweights`.


<a target='_blank' href='https://github.com/PierreBarrat/DCATools/blob/7b76e84e35cb66ec989505c4dced460cd09285ea/src/alignmenttools.jl#L41-L51' class='documenter-source'>source</a><br>


```
computefreqs(msa::String; q = 0, computew=false, weights=[], theta=0.2, saveweights="", header=false, format=1)
```

Compute pairwise frequencies for a file input. `msa` is a file containing the alignment in numerical format. Return frequencies `f1` and `f2` and weights.

Keywords: 

  * `q`: default `findmax(Y)[1]`
  * `weights`: default `[]`. If it is a `String`, phylogenetic weights are read from the corresponding file. If it is an `Array{Float64,1}`, they are used directly.
  * `computew`: default `false`. If true, phylogenetic weights are computed, calling the appropriate `computeweights`. `weights` is then ignored.
  * `saveweights` and `theta`: see `computeweights`.
  * `header` and `format`: see `readmsanum`.


<a target='_blank' href='https://github.com/PierreBarrat/DCATools/blob/7b76e84e35cb66ec989505c4dced460cd09285ea/src/alignmenttools.jl#L84-L95' class='documenter-source'>source</a><br>


<a id='Compute-weights-of-sequences.-1'></a>

### Compute weights of sequences.

<a id='DCATools.computeweights' href='#DCATools.computeweights'>#</a>
**`DCATools.computeweights`** &mdash; *Function*.



```
computeweights(msa::String; theta::Float64 = 0.2, saveweights::String = "", format=1, header=false)
```

Compute weights for file input. Compute weights of each sequence in file `msa` using reweighting threshold `theta` (default 0.2). 

Keywords: 

  * `theta`: threshold of similarity under which sequences are weighted down. Default `0.2`.
  * `saveweights`: weights are saved there if non empty. Default `""`
  * `format` and `header`: used to read `msa`. See `readmsanum`.


<a target='_blank' href='https://github.com/PierreBarrat/DCATools/blob/7b76e84e35cb66ec989505c4dced460cd09285ea/src/alignmenttools.jl#L106-L115' class='documenter-source'>source</a><br>


```
computeweights(Y::Array{Int64,2} ; theta = 0.2, saveweights="")
```

Compute weights for array input. Compute weights of each sequence in file `msa` using reweighting threshold `theta` (default 0.2). 

Keywords: 

  * `theta`: threshold of similarity under which sequences are weighted down. Default `0.2`.
  * `saveweights`: weights are saved there if non empty. Default `""`


<a target='_blank' href='https://github.com/PierreBarrat/DCATools/blob/7b76e84e35cb66ec989505c4dced460cd09285ea/src/alignmenttools.jl#L122-L130' class='documenter-source'>source</a><br>


```
computeweights(Y::Array{Int64,2},theta::Float64)
```

Basic routine. Compute weights of sequences in alignment `Y`, using threshold `theta`. 


<a target='_blank' href='https://github.com/PierreBarrat/DCATools/blob/7b76e84e35cb66ec989505c4dced460cd09285ea/src/alignmenttools.jl#L139-L143' class='documenter-source'>source</a><br>


<a id='Input-/-output-inputoutput.jl-1'></a>

## Input / output - inputoutput.jl


<a id='Read-DCA-parameters-from-file.-1'></a>

### Read DCA parameters from file.

<a id='DCATools.readparam' href='#DCATools.readparam'>#</a>
**`DCATools.readparam`** &mdash; *Function*.



```
readparam(infile::String ; format="mat", q=0)
```

Read dca parameters from `infile`. Format option can be either 

  * `"mcmc"`: `J i j a b val`
  * `"mat"`: One line of `infile` represents the vector `J[i,a][:]`. This is useful for parameters stored in dlm format. Optional argument `q` is needed in this case.

Output is of type `DCAgraph`.


<a target='_blank' href='https://github.com/PierreBarrat/DCATools/blob/7b76e84e35cb66ec989505c4dced460cd09285ea/src/inputoutput.jl#L7-L14' class='documenter-source'>source</a><br>


<a id='Write-DCA-parameters-to-file.-1'></a>

### Write DCA parameters to file.

<a id='DCATools.writeparam' href='#DCATools.writeparam'>#</a>
**`DCATools.writeparam`** &mdash; *Function*.



```
writeparam(outfile::String, g::DCAgraph; format="mat")
```

Write graph `g` to file `outfile`: 

  * as a matrix if `format=="mat"`
  * using `J i j a b value` if `format=="mcmc"`


<a target='_blank' href='https://github.com/PierreBarrat/DCATools/blob/7b76e84e35cb66ec989505c4dced460cd09285ea/src/inputoutput.jl#L98-L104' class='documenter-source'>source</a><br>


<a id='Read-a-numerical-alignment-from-a-file.-1'></a>

### Read a numerical alignment from a file.

<a id='DCATools.readmsanum' href='#DCATools.readmsanum'>#</a>
**`DCATools.readmsanum`** &mdash; *Function*.



```
readmsanum(infile::String ; format=1, header=false)
```

Read an MSA stored in `infile` in a numerical format. 

If `format=1`, amino acids should be mapped from 1 to `q`. If `format=0`, they should be mapped from 0 to `q-1`. `header` argument allows for discarding the first line of `infile`. 


<a target='_blank' href='https://github.com/PierreBarrat/DCATools/blob/7b76e84e35cb66ec989505c4dced460cd09285ea/src/inputoutput.jl#L142-L149' class='documenter-source'>source</a><br>


<a id='Model-tools-–-modeltools.jl-1'></a>

## Model tools – modeltools.jl


<a id='Switch-gauge-of-DCA-parameters.-1'></a>

### Switch gauge of DCA parameters.

<a id='DCATools.switchgauge!' href='#DCATools.switchgauge!'>#</a>
**`DCATools.switchgauge!`** &mdash; *Function*.



```
switchgauge!(g::DCAgraph ; gauge="0sum")
```

Switch parameters in `g` to gauge `gauge`.  Implemented gauges:

1. 0 sum: "0sum"
2. Lattice gas: "LG", "lg", "latticegas".


<a target='_blank' href='https://github.com/PierreBarrat/DCATools/blob/7b76e84e35cb66ec989505c4dced460cd09285ea/src/modeltools.jl#L3-L10' class='documenter-source'>source</a><br>


<a id='Compute-energies-of-sequences-for-given-DCA-parameters.-1'></a>

### Compute energies of sequences for given DCA parameters.

<a id='DCATools.computeenergies' href='#DCATools.computeenergies'>#</a>
**`DCATools.computeenergies`** &mdash; *Function*.



```
computeenergies(g::DCAgraph, sample::Array{Int64,2})
```

Compute energies of all configurations in `sample` with graph `g`.


<a target='_blank' href='https://github.com/PierreBarrat/DCATools/blob/7b76e84e35cb66ec989505c4dced460cd09285ea/src/modeltools.jl#L75-L79' class='documenter-source'>source</a><br>


```
computeenergies(g::DCAgraph, sample::Array{Int64,1})
```

Compute energies of all configurations in `sample` with graph `g`.


<a target='_blank' href='https://github.com/PierreBarrat/DCATools/blob/7b76e84e35cb66ec989505c4dced460cd09285ea/src/modeltools.jl#L95-L99' class='documenter-source'>source</a><br>


<a id='Infer-a-profile-model-from-frequencies-or-an-alignment.-1'></a>

### Infer a profile model from frequencies or an alignment.

<a id='DCATools.inferprofile' href='#DCATools.inferprofile'>#</a>
**`DCATools.inferprofile`** &mdash; *Function*.



```
inferprofile(Y::Array{Int64,2}; q=findmax(Y)[1], pc = 1e-5, weights=[], save::String="")
```

Infer profile model from alignment `Y`. 

Keywords:

  * `q`: Default to maximum value in `Y`.

-`pc` and `save`: See `inferprofile`. 

  * `weights`: see `computefreqs`
  * 


<a target='_blank' href='https://github.com/PierreBarrat/DCATools/blob/7b76e84e35cb66ec989505c4dced460cd09285ea/src/modeltools.jl#L105-L115' class='documenter-source'>source</a><br>


```
inferprofile(f1::Array{Float64,1}, q::Int64; pc = 1e-5, save::String="")
```

Infer profile model from frequencies `f1`. 

Keywords: -`pc`: Pseudocount ratio. Defaults to `1e-5`.

  * save: File to save inferred profile. Format of save is `mat`, readable by `readdlm` or `readparam`.


<a target='_blank' href='https://github.com/PierreBarrat/DCATools/blob/7b76e84e35cb66ec989505c4dced460cd09285ea/src/modeltools.jl#L121-L129' class='documenter-source'>source</a><br>


<a id='Compute-pseudo-likelihood-of-data-1'></a>

### Compute pseudo-likelihood of data

<a id='DCATools.pseudolikelihood' href='#DCATools.pseudolikelihood'>#</a>
**`DCATools.pseudolikelihood`** &mdash; *Function*.



```
pseudolikelihood(Y::Array{Int64,2}, g::DCAgraph; weights=ones(size(Y,1)))
```

Compute the pseudo-likelihood of configurations (*ie* sequences) in `Y` according to parameters in `g`. 


<a target='_blank' href='https://github.com/PierreBarrat/DCATools/blob/7b76e84e35cb66ec989505c4dced460cd09285ea/src/modeltools.jl#L143-L147' class='documenter-source'>source</a><br>


<a id='Contact-prediction-–-contactprediction.jl-1'></a>

## Contact prediction – contactprediction.jl


<a id='Compute-Fapc-score-for-given-coupling-matrix.-1'></a>

### Compute Fapc score for given coupling matrix.

<a id='DCATools.Fapc' href='#DCATools.Fapc'>#</a>
**`DCATools.Fapc`** &mdash; *Function*.



```
Fapc(A::Array{Float64,2}, q::Int64 ; APC::Bool = true, gap::Bool = false, cols::Int64=3)
```

Compute Frobenius norm of `q x q` blocks in matrix `A`. 

Keywords:

  * `APC`: apply the famous APC correction. Default to `true`.
  * `gap`: Remove the state `1` from the Frobenius norm. Default to `false`.
  * `col`: Format of output. With `3`, output is an array with rows being `i j value`. With `1`, output is a vector containing `value` only. Default to `3`.


<a target='_blank' href='https://github.com/PierreBarrat/DCATools/blob/7b76e84e35cb66ec989505c4dced460cd09285ea/src/contactprediction.jl#L47-L57' class='documenter-source'>source</a><br>


<a id='Compute-PPV-for-given-scores-and-distances.-1'></a>

### Compute PPV for given scores and distances.

<a id='DCATools.PPV' href='#DCATools.PPV'>#</a>
**`DCATools.PPV`** &mdash; *Function*.



```
PPV(scores::Array{Float64,2}, distances::Array{Float64,2} ; minrange=4, threshold=8)
```

Compute Positive Predictive Value (PPV) for `scores` and `distances`. Both of these arrays should be in the format `i j val`. 

Keywords:

  * `minrange`: Minimum value of `|i-j|` for a prediction to be made. Default to `4`.
  * `threshold`: Threshold defining contact. Default to `8`.


<a target='_blank' href='https://github.com/PierreBarrat/DCATools/blob/7b76e84e35cb66ec989505c4dced460cd09285ea/src/contactprediction.jl#L4-L12' class='documenter-source'>source</a><br>


<a id='Miscellaneous-–-misc.jl-1'></a>

## Miscellaneous – misc.jl


<a id='Fitting-quality-of-DCA-model-1'></a>

### Fitting quality of DCA model

<a id='DCATools.fitquality' href='#DCATools.fitquality'>#</a>
**`DCATools.fitquality`** &mdash; *Function*.



```
fitquality(f2_1::Array{Float64,2}, f1_1::Array{Float64,1}, f2_2::Array{Float64,2}, f1_2::Array{Float64,1}, q::Int64; withdiag=false)
```

Fitting quality of pairwise frequencies. Compare frequencies `(f1_1,f2_1)` (*e.g.* from natural sequences) to `(f1_2,f2_2)` (*e.g.* from a dca model). Output is 

1. Pearson correlation between connected correlations. (default: without diagonal elements)
2. Slope corresponding to the Pearson correlation.
3. Frobenius norm of the difference between connected correlations.
4. Same as 1. for magnetizations.
5. Same as 3. for magnetizations.


<a target='_blank' href='https://github.com/PierreBarrat/DCATools/blob/7b76e84e35cb66ec989505c4dced460cd09285ea/src/misc.jl#L3-L12' class='documenter-source'>source</a><br>


<a id='Three-points-correlation-1'></a>

### Three points correlation

<a id='DCATools.corr3p' href='#DCATools.corr3p'>#</a>
**`DCATools.corr3p`** &mdash; *Type*.



```
mutable struct corr3p
```

Three body correlation at positions `i, j, k` for states `a, b, c`, with value `cor`. 


<a target='_blank' href='https://github.com/PierreBarrat/DCATools/blob/7b76e84e35cb66ec989505c4dced460cd09285ea/src/misc.jl#L50-L54' class='documenter-source'>source</a><br>

<a id='DCATools.threepointscor' href='#DCATools.threepointscor'>#</a>
**`DCATools.threepointscor`** &mdash; *Function*.



```
threepointscor(Y::Array{Int64,2}, w::Array{Float64,1}, q::Int64, threshold::Float64)
```

Base routine. Compute three body correlations between columns of alignment `Y`, with sequences being weighted by `w`. Keep only correlation values that are above `threshold`. Return corresponding list of triplets `(i,j,k),(a,b,c)`. 


<a target='_blank' href='https://github.com/PierreBarrat/DCATools/blob/7b76e84e35cb66ec989505c4dced460cd09285ea/src/misc.jl#L65-L69' class='documenter-source'>source</a><br>


```
threepointscor(Y::Array{Int64,2}, w::Array{Float64,1}, q::Int64, triplets::Array{Int64,2})
```

Base routine. Compute three body correlations between columns of alignment `Y`, with sequences being weighted by `w`. Consider only triplets (i,j,k,a,b,c) in `triplets`. 

*Note*: this function could easily be optimized.


<a target='_blank' href='https://github.com/PierreBarrat/DCATools/blob/7b76e84e35cb66ec989505c4dced460cd09285ea/src/misc.jl#L116-L122' class='documenter-source'>source</a><br>


```
threepointscor(Y::Array{Int64,2}; q = findmax(Y)[1], threshold=0, triplets = Array{Int64,2}(undef,0,1),computew=false, weights=[], theta=0.2, saveweights="")
```

Compute three points frequencies for an array input. `Y` is an array of `Int64`. Return an array of `c3p` structures. 

Keywords: 

  * `q`: default `findmax(Y)[1]`
  * `threshold`: default `0`. Only points with correlation higher than `threshold` are returned.
  * `triplets`: default `Array{Int64,2}(undef,0,1)`. List of position and amino acid triplets `(i,j,k,a,b,c)` for which the correlation is computed.
  * `weights`: default `[]`. If it is a `String`, phylogenetic weights are read from the corresponding file. If it is an `Array{Float64,1}`, they are used directly.
  * `computew`: default `false`. If true, phylogenetic weights are computed, calling the appropriate `computeweights`. `weights` is then ignored.
  * `saveweights` and `theta`: see `computeweights`.


<a target='_blank' href='https://github.com/PierreBarrat/DCATools/blob/7b76e84e35cb66ec989505c4dced460cd09285ea/src/misc.jl#L155-L167' class='documenter-source'>source</a><br>

