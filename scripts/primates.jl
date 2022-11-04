using RCall
using Distributions
using StatsPlots
using Diversification

##############################
##
##   Load the data files
##
###############################
treefile = "data/primates.tre"
datafile = ""
phy = readtree(treefile)
num_total_species = 367
ρ = length(phy[:tip_label]) / num_total_species
data = make_SSEdata(phy, datafile, ρ; include_traits = false)

##############################
##
##   Set up the model
##
###############################
λ = [0.05, 0.1, 0.15, 0.2]
μ = 0.5 .* λ

η = 0.1
model = SSEconstant(λ, μ, η)


## Calculate the backwards-forwards pass equations
Ds, Fs = backwards_forwards_pass(model, data; verbose = true) 
res = calculate_tree_rates(data, model, Ds, Fs; verbose = false);

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
    geom_tiplab(size = 8) +
    theme(legend.position = c(0.2, 0.8)) +
    xlim(c(0.0, th + 10)) 
""";

