using RCall
using Distributions
using StatsPlots
using Diversification
using Optim

##############################
##
##   Load the data files
##
###############################
treefile = "data/primates.tre"
phy = readtree(treefile)
num_total_species = 367
ρ = length(phy[:tip_label]) / num_total_species
data = make_SSEdata2(phy, ρ)

##############################
##
##   Set up the model
##
###############################

H = 0.587405

d1 = LogNormal(log(0.29), 3 * H)
d2 = LogNormal(log(0.20), 3 * H)
η = 0.1

n = 10
speciation = make_quantiles(d1, n)
extinction = make_quantiles(d2, n)

k = n ^2
λ = zeros(k)
μ = zeros(k)

for (i, (sp, ex)) in enumerate(Iterators.product(speciation, extinction))
    λ[i] = sp
    μ[i] = ex
end
model = SSEconstant(λ, μ, η)

##############################
##
## Calculate the backwards-forwards pass equations
##
###############################


Ds, Fs = backwards_forwards_pass(model, data; verbose = true);
Ps = ancestral_state_probabilities(data, model, Ds, Fs);
res = calculate_tree_rates(data, model, Ds, Fs, Ps; verbose = false);

average_node_rates = res["average_node_rates"]

phy = Dict("edge" => data.edges,
      "tip.label" => data.tiplab,
      "Nnode" => length(data.tiplab)-1,
     "edge.length" => data.branch_lengths)

lambda_average = average_node_rates["λ"]

##############################
##
## Plot the tree using some R code
##
###############################

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

