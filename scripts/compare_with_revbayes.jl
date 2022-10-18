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
datafile = "data/bears.csv"
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
Ds, Fs = backwards_forwards_pass(model, data; verbose = false) 
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

df1 <- tibble("node" = 1:max(phy$edge),
            "Speciation rate" = lambda_average)
x <- as_tibble(phy)

phydf <- merge(x, df1, by = "node")
td_phy <- as.treedata(phydf)

p1a <- ggtree(td_phy, aes(color = `Speciation rate`)) +
  ggtitle("backwards-forwards approach")
p1 <- p1a +
  geom_tiplab()
ggsave("figures/p1a.pdf", p1a, width = 800, height = 800, units = "mm")
"""

## Compare with RevBayes output
R"""
source("scripts/matchnodes.R")

mn <- matchNodes(phy)

rates <- as_tibble(read.table("output/bears_BDS_rates.log", header = TRUE))
rates <- select(rates, starts_with("avg_lambda"))

lambda_means <- unname(sapply(rates, mean))
my_order <- mn[order(mn$Rev),"R"]
lambda_means <- lambda_means[my_order]

d2 <- tibble("node" = 1:max(phy$edge),
            "Speciation rate" = lambda_means)
phydf2 <- merge(x, d2, by = "node")
phy2 <- as.treedata(phydf2)
p2 <- ggtree(phy2, aes(color = `Speciation rate`)) +
    geom_tiplab() +
    ggtitle("RevBayes")

p1 | p2 
"""

@rget rates
@rget mn
mn = convert.(Int64, mn)

vplot = StatsPlots.plot(ylab = "Average speciation rate per branch", xlab = "branch index (ape node index)")
for (i, Rev_index) in enumerate(mn[!,"Rev"])
    if i == 1
        l1 = "RevBayes λ"
        l2 = "RevBayes median"
        l3 = "Revbayes mean"
    else
        l1 = ""
        l2 = ""
        l3 = ""
    end
    StatsPlots.violin!(vplot, [i], [rates[!, Rev_index]], label = l1, color = "lightblue")
#    StatsPlots.scatter!(vplot, [i], [median(rates[!, Rev_index])], color = "orange", label = l2, alpha = 0.5)
    StatsPlots.scatter!(vplot, [i], [mean(rates[!, Rev_index])], color = "red", label = l3, alpha = 0.7)
end
StatsPlots.scatter!(vplot, 1:15, average_node_rates["λ"], label = "New, Backwards-forwards pass", color = "black", alpha = 0.7)
vplot2 = deepcopy(vplot)
plot!(vplot2, ylim = (0.0, 0.5), title = "different y-axis limits")
branch_rateplot = plot(vplot, vplot2)

savefig(branch_rateplot, "figures/branch_rates_RevBayes_versus_new_approach_bears.pdf")
