#bamboo

using Distributions
using Diversification

phy = readtree("data/152_plastid_RAxML_bipartitions.tre")
ntotal = 1450
ρ = length(phy[:tip_label]) / ntotal
data = make_SSEdata2(phy, ρ)

