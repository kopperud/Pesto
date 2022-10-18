using CSV
using RCall
using DifferentialEquations
using Distributions
using ProgressMeter
using StatsPlots

include("datatypes.jl")
include("utils.jl")
include("diversitree_bisse.jl")
include("postorder.jl")
include("preorder.jl")
include("ODE.jl")
include("logLroot.jl")
include("constantBDP.jl")
include("tree_rates.jl")
include("backwards_forwards.jl")

##############################
##
##   Load the data files
##
###############################
treefile = "data/bears.tre"
datafile = "data/bears.csv"
#treefile = "data/actinopterygii/actinopt_full_1.tree"
#datafile = ""
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

bs = []
ks = [2, 4, 8, 16, 32, 64, 128, 256]
for k in ks 
    λ = make_quantiles(d1, k)
    μ = make_quantiles(d2, k)
    η = 0.1
    k = length(λ)

    model = SSEconstant(λ, μ, η)

    ## Calculate the backwards-forwards pass equations
    b = @benchmark Ds, Fs = backwards_forwards_pass(model, data; verbose = false) samples = 50
    Ds, Fs = backwards_forwards_pass(model, data; verbose = false) 
    res = calculate_tree_rates(data, model, Ds, Fs; verbose = false);
    average_node_rates = res["average_node_rates"]

    phy = Dict("edge" => data.edges,
              "tip.label" => data.tiplab,
              "Nnode" => length(data.tiplab)-1,
             "edge.length" => data.branch_lengths)
    append!(bs, [b])
    print(".")
end

tplot = plot()
for (k, b) in zip(ks, bs)
    violin!(tplot, [k], [b.times .* 10^(-9)], lab = string(k)*" categories")
end
plot!(tplot, [ks[1], ks[end]], [mean(bs[1].times), mean(bs[end].times)] .* 10^(-9), lab = "Linear scaling")
plot!(tplot, legend = :topleft, ylab = "seconds (50 replicates)")


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
#ggsave("figures/p1.pdf", p1a, width = 800, height = 800, units = "mm")
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
plot(vplot, vplot2)



ps = []
for i in 1:14
    times1 = range(minimum(Fs[i].t), maximum(Fs[i].t); length = 100)
    p1 = plot(times1, hcat(Ds[i].(times1)...)', title = "D")
    p2 = plot(times1, hcat(Fs[i].(times1)...)', title = "F")
    p3 = plot(times1, hcat(Ps[i].(times1)...)', title = "P")
    p = plot(p1, p2, p3, layout = (1,3), xflip = true, nrow = 1 )
    append!(ps, [p])
end


#if false
#    println("Diversitree-BISSE model:\n")
#    asr = anc_state_prob_bisse(treefile, datafile, model)
#    println("Ancestral state marginal probabilities:")
#    display(asr)
#end


