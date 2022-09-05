using Test
using DCATools

mkpath("../tmp/")

L = 5
q = 4
g = DCAGraph(L, q; init = :rand)

t = @testset "Extended format - index 1" begin
	write("../tmp/gext.dat", g, :extended; index_style=1, sigdigits=5)
	gg = DCAGraph("../tmp/gext.dat")
	@test g.L == gg.L
	@test g.q == gg.q
	@test isapprox(g.J, gg.J, rtol=1e-4)
	@test isapprox(g.h, gg.h, rtol=1e-4)
end

t = @testset "Extended format - index 0" begin
	write("../tmp/gext.dat", g, :extended; index_style=0, sigdigits=5)
	gg = DCAGraph("../tmp/gext.dat")
	@test g.L == gg.L
	@test g.q == gg.q
	@test isapprox(g.J, gg.J, rtol=1e-4)
	@test isapprox(g.h, gg.h, rtol=1e-4)
end

t = @testset "Matrix" begin
	write("../tmp/gmat.dat", g, :matrix; sigdigits=5)
	gg = DCAGraph("../tmp/gmat.dat", :matrix; q)
	@test g.L == gg.L
	@test g.q == gg.q
	@test isapprox(g.J, gg.J, rtol=1e-4)
	@test isapprox(g.h, gg.h, rtol=1e-4)
end



