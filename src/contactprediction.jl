export PPV, Fapc


"""
    PPV(scores::Array{Float64,2}, distances::Array{Float64,2} ; minrange=4, threshold=8)

Compute Positive Predictive Value (PPV) for `scores` and `distances`. Both of these arrays should be in the format `i j val`. 

Keywords:
- `minrange`: Minimum value of `|i-j|` for a prediction to be made. Default to `4`.
- `threshold`: Threshold defining contact. Default to `8`. 
"""
function PPV(scores::Array{Float64,2}, distances::Array{Float64,2} ; minrange=4, threshold=8)
    if size(scores,1)!=size(distances,1)
        error("contactprediction.jl - PPV: `scores` and `distances` do not have the same size.")
    elseif size(scores[1,:],1)!=3 || size(scores[1,:],1)!=3
        error("contactprediction.jl - PPV: `scores` and `distances` should be in the format `i j value`.")
    end

    scores = sortslices(scores, dims=1, rev=true, by=x->x[3])
    d = Dict{Tuple{Float64,Float64},Float64}()
    for l in 1:size(distances,1)
        d[(distances[l,1], distances[l,2])] = distances[l,3]
    end
    tp = 0
    np = 0
    TP = []
    cc=0
    for r in 1:size(scores,1)
        if abs(scores[r,1] - scores[r,2]) > minrange && get(d, (scores[r,1], scores[r,2]), -1.)>=0.
            np +=1
            if d[(scores[r,1], scores[r,2])]<threshold 
                tp+=1
            else
                # cc<10?println("$(scores[r,1]), $(scores[r,2])"):Void
                cc+=1
            end
            push!(TP,tp/np)
        end
    end

    return TP
end



"""
    Fapc(A::Array{Float64,2}, q::Int64 ; APC::Bool = true, gap::Bool = false, cols::Int64=3)

Compute Frobenius norm of `q x q` blocks in matrix `A`. 

Keywords:
- `APC`: apply the famous APC correction. Default to `true`.
- `gap`: Remove the state `1` from the Frobenius norm. Default to `false`.
- `col`: Format of output. With `3`, output is an array with rows being `i j value`. With `1`, output is a vector containing `value` only. Default to `3`. 

"""
function Fapc(A::Array{Float64,2}, q::Int64 ; APC::Bool = true, gap::Bool = false, cols::Int64=3)

    (L,L) = size(A)
    try
        L = convert(Int64,L/q)
    catch
        println("Size of input matrix is ",size(A))
        error("contactprediction.jl - Fapc: Input matrix A does not have the correct size.")
    end
    if !(A==A')
        error("contactprediction.jl - Fapc: Input matrix A is not symmetric")
    end

    S = zeros(Float64, L,L)
    for i in 0:L-1
        for j in (i+1):L-1
            S[i+1,j+1] = gap ? sqrt(sum(A[i*q .+ (1:q), j*q .+ (1:q)].^2)) : sqrt(sum(A[i*q .+ (2:q), j*q .+ (2:q)].^2))
            S[j+1,i+1] = S[i+1,j+1]
        end
    end

    if cols == 1
        F = zeros(Float64, convert(Int64, L*(L-1)/2))
        pos = 1
        for i in 1:L
            for j in (i+1):L
                F[pos] = S[i,j] - (mean(S[:,j]) * mean(S[:,i]) / mean(S)) * APC
                pos += 1
            end
        end
    elseif cols ==3 
        F = zeros(Float64, convert(Int64, L*(L-1)/2),3)
        pos = 1
        for i in 1:L
            for j in (i+1):L
                F[pos,:] = [i,j,S[i,j] - (mean(S[:,j]) * mean(S[:,i]) / mean(S)) * APC]
                pos += 1
            end
        end
    end
    return F
end