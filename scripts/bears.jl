using RCall
using Distributions
using StatsPlots
using Diversification

##############################
##
##   Load the data files
##
###############################
treefile = "data/bears.tre"
datafile = ""
phy = readtree(treefile)
num_total_species = 8
ρ = length(phy[:tip_label]) / num_total_species
data = make_SSEdata(phy, datafile, ρ; include_traits = false)

##############################
##
##   Set up the model
##
###############################
λ = [0.1, 0.2]
μ = [0.05, 0.15]

η = 0.05
model = SSEconstant(λ, μ, η)


## Calculate the backwards-forwards pass equations
Ds, Fs = backwards_forwards_pass(model, data; verbose = true) 
Ps = ancestral_state_probabilities(data, model, Ds, Fs)
res = calculate_tree_rates(data, model, Ds, Fs, Ps; verbose = false);

average_node_rates = res["average_node_rates"]

phy = Dict("edge" => data.edges,
      "tip.label" => data.tiplab,
      "Nnode" => length(data.tiplab)-1,
     "edge.length" => data.branch_lengths)

lambda_average = average_node_rates["λ"]

@rput lambda_average
@rput phy
R"""

class(phy) <- "phylo"
th <- max(ape::node.depth.edgelength(phy))

df1 <- tibble::tibble("node" = 1:max(phy$edge),
            "Speciation rate" = lambda_average)
x <- tidytree::as_tibble(phy)

phydf <- merge(x, df1, by = "node")
td_phy <- tidytree::as.treedata(phydf)

p1a <- ggtree::ggtree(td_phy, aes(color = `Speciation rate`)) +
    ggtree::geom_tiplab(size = 8) +
    ggplot2::theme(legend.position = c(0.2, 0.8)) +
    ggplot2::xlim(c(0.0, th + 10)) 
""";

