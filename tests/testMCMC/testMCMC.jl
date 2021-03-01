using Test
using DCATools, DCATools.MCMC
using Statistics, Polynomials

q = 21
g = readparam("testMCMC/parameters_PF01535_mat.txt", q=q)
sample_ref = readmsanum("testMCMC/MC_matteo.txt", format=0, header=true)
f1,f2 = computefreqs(sample_ref)

sample_test = doMCMC(g, 5_000, 5)
p1,p2 = computefreqs(sample_ref)

@testset "Correlation/slope of frequencies" begin
	P1 = fit(f1,p1,1)
	P2 = fit(vec(f2), vec(p2), 1)
	@test isapprox(cor(f1,p1), 1., atol=1e-3)
	@test isapprox(cor(vec(f2),vec(p2)), 1., atol=1e-3)
	@test isapprox(P1[1], 1., atol=1e-3)
	@test isapprox(P2[1], 1., atol=1e-3)
end
