using DCATools, DCATools.BM

aln = readmsanum("PF00096/msa_PF00096.txt", format=0, header=true)[1:5000,:];
f1, f2, w = pairwise_frequencies(aln,computew=true);
L = size(aln)[2]
q = 21

ginit = inferprofile(f1, q)
out = bmlearn(f1,f2,L,q, ginit=ginit, savefolder="testingBM", verbose=true, nit=10, comment="Simple test run")
