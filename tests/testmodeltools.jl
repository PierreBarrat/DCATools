using DCATools
using Test
using Polynomials, Statistics

N = 10 
q = 4 
M = 1000

J = rand(Float64, N*q,N*q)
J = J + J'
for i = 1:N
    J[(i-1)*q .+ (1:q), (i-1)*q .+ (1:q)] .= 0
end
h = rand(Float64,N*q)

g = DCAgraph(J,h,N,q)
g0 = deepcopy(g)
glg = deepcopy(g)
switchgauge!(g0)
switchgauge!(glg, gauge=:lg)

@testset "0-sum gauge" begin
	@test isapprox(sum(g0.J), 0., atol = 1e-12)
	@test isapprox(sum(g0.h), 0., atol = 1e-12)
	@test isapprox(findmax(sum(g0.J))[1], 0., atol=1e-12)
	@test isapprox(sum(g0.h), 0., atol=1e-12)
end

@testset "Lattice gas gauge" begin
	@test isapprox(findmax(glg.J[q,:])[1], 0., atol=1e-12)
	@test isapprox(findmax(glg.J[:,q])[1], 0., atol=1e-12)
	@test isapprox(findmax(glg.h[q:q:end])[1], 0., atol=1e-12)
end 

@testset "Gauge change: conservation of energies" begin
	seq = rand(collect(1:q), M, N)
	E_init = computeenergies(g,seq)
	E_0 = computeenergies(g0,seq)
	E_lg = computeenergies(glg, seq)
	P1 = fit(E_init, E_0, 1)
	P2 = fit(E_init, E_lg, 1)
	@test isapprox(P1[1], 1, atol=1e-3)
	@test isapprox(P2[1], 1, atol=1e-3)
	@test isapprox(cor(E_init, E_0), 1, atol=1e-3)
	@test isapprox(cor(E_init, E_lg), 1, atol=1e-3)
end
