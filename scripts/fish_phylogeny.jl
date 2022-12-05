using RCall
using Distributions
using StatsPlots
using Diversification

##############################
##
##   Load the data files
##
###############################
treefile = "data/actinopterygii/actinopt_full_1.tree"
datafile = ""
phy = readtree(treefile)
ρ = 1.0
data = make_SSEdata(phy, datafile, ρ; include_traits = false)

##############################
##
##   Set up the model
##
###############################
H = 0.587405
d1 = LogNormal(log(0.08251888), 3 * H)
d2 = LogNormal(log(0.08251888/2.0), 3 * H)

k = 20
λ = make_quantiles(d1, k)
μ = make_quantiles(d2, k)
η = 0.1

model = SSEconstant(λ, μ, η)

## Calculate the backwards-forwards pass equations
Ds, Fs = backwards_forwards_pass(model, data; verbose = true) 
Ps = ancestral_state_probabilities(data, model, Ds, Fs)
res = calculate_tree_rates(data, model, Ds, Fs, Ps; verbose = true);
average_node_rates = res["average_node_rates"]

phy = Dict("edge" => data.edges,
      "tip.label" => data.tiplab,
      "Nnode" => length(data.tiplab)-1,
     "edge.length" => data.branch_lengths)

lambda_average = average_node_rates["λ"]
@rput lambda_average
@rput phy
R"""
library(ape)
library(ggtree)
library(tidytree)
library(ggplot2)
library(dplyr)

class(phy) <- "phylo"
th <- max(node.depth.edgelength(phy))

df1 <- tibble("node" = 1:max(phy$edge),
            "Speciation rate" = lambda_average)
x <- as_tibble(phy)

phydf <- merge(x, df1, by = "node")
td_phy <- as.treedata(phydf)

p1a <- ggtree(td_phy, aes(color = `Speciation rate`)) +
  ggtitle("backwards-forwards approach")
p1 <- p1a +
  geom_tiplab() +
    theme(legend.position = c(0.2, 0.8)) +
    xlim(c(0.0, th + 10))

ggsave("figures/fish_rates.pdf", p1a, width = 800, height = 800, units = "mm")
"""

