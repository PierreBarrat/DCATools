print("\n\n --- Testing aligntools.jl ---\n\n")

using DCATools
using DelimitedFiles
using Test

## Testing compute freqs
#=
# Cases to test: 
# - With a numerical alignment. 
# 1. No weights.
# 2. Weights passed as array.
# 3. Weights passed as a file.
# 4. Weights to be computed, not saved.
# 5. Weights to be computed and saved. 
# - With the alignment file 
=#

fw_2 = readdlm("testaligntools/freqs_w.txt",Float64)
fw_1 = fw_2[end,:]
fw_2 = fw_2[1:(end-1),:]
fnw_2 = readdlm("testaligntools/freqs_nw.txt", Float64)
fnw_1 = fnw_2[end,:]
fnw_2 = fnw_2[1:(end-1),:]

ww = [1.,1.]
wnw = [0.5, 0.5, 1.]

Yu = readmsanum("testaligntools/msa_n3q2_unique.txt")
Ynu = readmsanum("testaligntools/msa_n3q2_nonunique.txt")

# 1.
@testset "No weights" begin
	f1,f2 = pairwise_frequencies(Yu)
	g1,g2 = pairwise_frequencies(Ynu)
	@test f1 == fw_1
	@test f2 == fw_2
	@test g1 == fnw_1
	g2 == fnw_2
end

# 2. 
@testset "Weights passed as array" begin
	f1,f2 = pairwise_frequencies(Yu,weights=ww)
	g1,g2 = pairwise_frequencies(Ynu,weights=wnw)
	@test f1==fw_1 && f2==fw_2 && g1==fw_1 && g2==fw_2
end

# 3. 
@testset "Weights passed as file" begin
	f1,f2 = pairwise_frequencies(Yu,weights="testaligntools/weights_w.txt")
	g1,g2 = pairwise_frequencies(Ynu,weights="testaligntools/weights_nw.txt")
   	@test f1==fw_1 && f2==fw_2 && g1==fw_1 && g2==fw_2
end

# 4. 
@testset "Weights to be computed and not saved" begin
	@test_logs (:warn,"both keywords `weights` and `computew` declared. `weights` ignored.\n") pairwise_frequencies(Yu,computew=true, weights="toto.txt")
	f1,f2,w = pairwise_frequencies(Yu,computew=true)
	g1,g2,nw = pairwise_frequencies(Ynu,computew=true)
	@test f1==fw_1 && f2==fw_2 && g1==fw_1 && g2==fw_2 && w==ww && nw==wnw
end

# 5. 
@testset "Weights to be computed and saved" begin
	f1,f2,w = pairwise_frequencies(Yu,computew=true, saveweights="testaligntools/testsave_ww.txt")
	g1,g2,nw = pairwise_frequencies(Ynu,computew=true, saveweights="testaligntools/testsave_nw.txt")
	@test f1==fw_1 && f2==fw_2 && g1==fw_1 && g2==fw_2 && w==ww && nw==wnw
end

# 6.
@testset "With alignment file 1." begin
	f1,f2 = pairwise_frequencies("testaligntools/msa_n3q2_unique.txt")
	g1,g2 = pairwise_frequencies("testaligntools/msa_n3q2_nonunique.txt")
	@test f1==fw_1 && f2==fw_2 && g1==fnw_1 && g2==fnw_2
end

# 7.
@testset "With alignment file 2." begin
	f1,f2, w = pairwise_frequencies("testaligntools/msa_n3q2_unique.txt",computew=true, saveweights="testaligntools/testsave_ww.txt")
	g1,g2, nw = pairwise_frequencies("testaligntools/msa_n3q2_nonunique.txt",computew=true, saveweights="testaligntools/testsave_nw.txt")
	@test f1==fw_1 && f2==fw_2 && g1==fw_1 && g2==fw_2 && w==ww && nw==wnw
end
