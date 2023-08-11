"""
	mf_couplings(C ; q = 21)

never tested, I would be surprised if it works
"""
function mf_couplings(C ; q = 21)
	L = Int64(size(C)[1]/q)
	println(L)
	Cr = remove_column(C, L, q, q)
	iC = inv(Cr)
	J = zeros(L*q, L*q)
	for i in 1:L
		J[(i-1)*q.+(1:q), (i-1)*q.+(1:q)] .= 0
		for j in (i+1):L
			J[(i-1)*q.+(1:q-1), (j-1)*q.+(1:q-1)] .= iC[(i-1)*(q-1) .+ (1:q-1), (j-1)*(q-1) .+ (1:q-1)]
		end
	end
	J .+= J'
	return -J
end

"""
	remove_column(A::Array{T,2}, L, q, col) where T
"""
function remove_column(A::Array{T,2}, L, q, col) where T
	Ar = zeros(T, L*(q-1), L*(q-1))
	for i in 1:L
		for j in (i+1):L
			Ar[(i-1)*(q-1) .+ (1:q-1), (j-1)*(q-1) .+ (1:q-1)] = A[(i-1)*q .+ (1:q-1), (j-1)*q .+ (1:q-1)]
			Ar[(j-1)*(q-1) .+ (1:q-1), (i-1)*(q-1) .+ (1:q-1)] = Ar[(i-1)*(q-1) .+ (1:q-1), (j-1)*(q-1) .+ (1:q-1)]'
		end
		Ar[(i-1)*(q-1) .+ (1:q-1), (i-1)*(q-1) .+ (1:q-1)] = A[(i-1)*q .+ (1:q-1), (i-1)*q .+ (1:q-1)]
	end
	return Ar
end
