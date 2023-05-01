using Test
using DCATools
using Statistics, Polynomials

q = 21
g = DCAGraph("MCMC/parameters_PF01535_mat.txt", :matrix; q=q)
sample_ref = read_msa_num("MCMC/MC_matteo.txt", index_style=0, header=true)
f1, f2 = pairwise_frequencies(sample_ref)

sample_test = sample(g, 5_000; Twait=5)
p1, p2 = pairwise_frequencies(sample_ref)

@testset "Correlation/slope of frequencies" begin
	P1 = fit(f1, p1,1)
	P2 = fit(vec(f2), vec(p2), 1)
	@test isapprox(cor(f1,p1), 1., atol=1e-3)
	@test isapprox(cor(vec(f2),vec(p2)), 1., atol=1e-3)
	@test isapprox(P1[1], 1., atol=1e-3)
	@test isapprox(P2[1], 1., atol=1e-3)
end
