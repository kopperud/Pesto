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


ggsave("figures/p1a.pdf", p1a, width = 800, height = 800, units = "mm")
"""

## Compare with RevBayes output
R"""
source("scripts/matchnodes.R")

mn <- matchNodes(phy)

rates <- as_tibble(read.table("output/bears_BDS_rates.log", header = TRUE))
rates <- select(rates, starts_with("avg_lambda"))

lambda_means <- sapply(1:max(phy$edge), function(i) mean(rates[[i]]))
my_order <- mn[order(mn$Rev),"R"]

d2 <- tibble("node" = my_order,
            "Speciation rate" = lambda_means)
d2[d2$node == length(phy$tip.label)+1, "Speciation rate"] <- NA ## assign NA at root
phydf2 <- merge(x, d2, by = "node")
phy2 <- as.treedata(phydf2)
p2 <- ggtree(phy2, aes(color = `Speciation rate`)) +
    geom_tiplab() +
    ggtitle("RevBayes") +
    theme(legend.position = c(0.2, 0.8)) +
    xlim(c(0.0, th+10))

p_combined <- p2 | p1 
ggsave("figures/tree_rb_fb_pdf.pdf", p_combined, width = 250, height = 250, units = "mm")
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
    StatsPlots.scatter!(vplot, [i], [mean(rates[!, Rev_index])], color = "red", label = l3, alpha = 0.7)
end
StatsPlots.scatter!(vplot, 1:15, average_node_rates["λ"], label = "New, Backwards-forwards pass", color = "black", alpha = 0.7)
vplot2 = deepcopy(vplot)
plot!(vplot2, ylim = (0.0, 0.5), title = "different y-axis limits")
branch_rateplot = plot(vplot, vplot2)

savefig(branch_rateplot, "figures/branch_rates_RevBayes_versus_new_approach_bears.pdf")

Ps = res["Ps"]
Fs = res["Fs"]
Ds = res["Ds"]

ps = []
for i in 1:14
    times1 = range(minimum(Fs[i].t), maximum(Fs[i].t); length = 100)
    p1 = plot(times1, hcat(Ds[i].(times1)...)', title = "D")
    p2 = plot(times1, hcat(Fs[i].(times1)...)', title = "F")
    p3 = plot(times1, hcat(Ps[i].(times1)...)', title = "P")
    p = plot(p1, p2, p3, layout = (1,3), xflip = true, nrow = 1 )
    append!(ps, [p])
end

## histogram of rates
hs = []
for (i, col) in enumerate(eachcol(rates))
    p = histogram(col, xlab = "rate", ylab = "frequency", title = "branch "*string(i))
    append!(hs, [p])
end
savefig(plot(hs[1:4]...), "figures/stochastic_character_map_sample_bears.pdf")

xs = []
ys = average_node_rates["λ"]
xerrors1 = []
xerrors2 = []
for (i, Rev_index) in enumerate(mn[!,"Rev"])
    x = mean(rates[!, Rev_index])
    xerror = quantile(rates[!, Rev_index], [0.025, 0.975])

    append!(xs, x)
    append!(xerrors1, xerror[1])
    append!(xerrors2, xerror[2])
end

cplot = scatter(xs, ys, xerror = (xerrors1, xerrors2), xlab = "RevBayes (mean + 95% CI)", ylab = "New approach", lab = "Average branch speciation rate", legend = :bottomright)
plot!(cplot, [0.05, 0.25], [0.05, 0.25], color = "gray", lab = "One-to-one")


## varying time points
using DataFrames
using CSV

dfs = Dict()
nts = [500, 1000, 1500, 5000, 10000]
for nt in nts
    df = CSV.read("output/bears_BDS_rates_" * string(nt)* "_run_1.log", DataFrame)
    dfs[nt] = df
end

l = Dict()
for nt in nts
    df = dfs[nt]

    res = []
    for (i, Rev_index) in enumerate(mn[!, "Rev"])
        append!(res, mean(df[!, "avg_lambda[" * string(Rev_index) * "]"]))
    end
    l[nt] = res
end

y = hcat(values(l)...)

ps2 = []
for i in 1:15
    p = plot(nts, y[i,:], title = "branch "*string(i), lab = "", xlab = "Number of time slices", ylab = "Mean rate", xscale = :log)
    scatter!(p, nts, y[i,:], lab = "RevBayes")
    plot!(p, [nts[1], nts[end]], [average_node_rates["λ"][i], average_node_rates["λ"][i]], lab = "Julia impl")
    append!(ps2, [p])
end


p_timeslices = plot(ps2[1:6]..., size = (800, 800))
savefig(p_timeslices, "figures/p_timeslices.pdf")






