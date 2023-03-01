using Diversification
using Distributions
using Plots
using ProgressMeter
using BenchmarkTools
using RCall
using DataFrames
using CSV
using StatsBase
#using LaTeXStrings
#using Measures
using LinearAlgebra

phy = readtree("data/cetaceans.tre")
ρ = 1.0
data = make_SSEdata2(phy, ρ)

#λ = [0.1, 0.3, 0.1, 0.3]
#μ = [0.05, 0.10, 0.15, 0.20]
λ = [0.2, 0.1]
μ = [0.1, 0.0]
η = 0.1
K = length(λ)

model = SSEconstant(λ, μ, η)

res = birth_death_shift(model, data)

Ds, Fs = Diversification.backwards_forwards_pass(model, data);
Ss = Diversification.ancestral_state_probabilities(data, model, Ds, Fs)
E = extinction_probability(model, data)

## revbayes to ape indices
R"""
library(ape)
phy <- read.tree("data/cetaceans.tre")
source("scripts/matchnodes.R")
mn <- matchNodes(phy)
"""
@rget mn
mn[!, :Rev] = convert(Vector{Int64}, mn[!, :Rev])
Rev_to_R = Dict()
R_to_Rev = Dict()
for row in eachrow(mn)
    Rev_to_R[row["Rev"]] = row["R"]
    R_to_Rev[row["R"]] = row["Rev"]
end

## load the log file
## load the RevBayes log files
#df = DataFrame(CSV.File("output/BDS_SCM_cetaceans_2.log"))
df = DataFrame(CSV.File("output/cetaceans_scm_2500.log"))

nbranches = length(data.branch_lengths)
Nscm = zeros(nbranches)
root_index = length(data.tiplab) + 1
for edge_idx in 1:nbranches
    node_index = data.edges[edge_idx, 2]

    if node_index != root_index
        Rev_index = R_to_Rev[node_index]
        Nscm[edge_idx] = mean(df[!, "num_shifts["*string(Rev_index)*"]"])
    end
end

#ntimes = [10, 10, 25, 50, 75, 100, 150, 200, 250, 500, 1000, 2500]

ntimes = Int64.(round.(collect(range(5, 2500; length=10))))

ntimes = [50, 100, 150, 250, 350, 500, 1000]
nshiftsx = zeros(length(ntimes))
@showprogress for (i, nt) in enumerate(ntimes)
    nshiftsx[i] = sum(compute_nshifts(model, data, Ds, Ss; ntimes=nt))
end



p = plot(ntimes, nshiftsx, seriestype=[:scatter],
    color="black", xlab="Number of time slices\n(per branch)",
    ylab="Sum of E[N] across all branches", label="Analytical solution\n(difference equation)", grid=false)
plot!(p, ntimes, nshiftsx, seriestype=:line, color="black", label="")
plot!(p, [ntimes[1], ntimes[end]], [sum(Nscm), sum(Nscm)], color="orange",
    label="Stochastic character map\n(2000 iters)")




nshifts = compute_nshifts(model, data, Ds, Ss; ntimes=1000, ape_order=false)

sum(nshifts)


histogram(nshifts, bins=15)
std_err_mean = sqrt.(Nscm ./ size(df)[1])
shiftplot = plot(nshifts, Nscm, yerror=std_err_mean, linetype=:scatter,
    xlab="Analytical solution", grid=false,
    ylab="RevBayes SCM", lab="Number of shifts")
ymax = maximum([maximum(Nscm), maximum(nshifts)])
plot!(shiftplot, [0.0, ymax], [0.0, ymax],
    lab="One-to-one", linestyle=:dash)
savefig(shiftplot, "figures/scm_vs_nshift.pdf")

##

plot(nshifts, data.branch_lengths .* η, linetype=:scatter, xlab="nshifts", ylab="bl * η")
plot!([0.0, 0.2], [0.0, 0.3])
plot!([nshifts[104]], [data.branch_lengths[104]], color="red", lab="", linetype=:scatter)
annotate!(nshifts, data.branch_lengths, 1:nbranches)

## reorder to ape indices
node_nshifts = compute_nshifts(model, data, Ds, Ss; ntimes=80)

phy = Dict("edge" => data.edges,
    "tip.label" => data.tiplab,
    "Nnode" => length(data.tiplab) - 1,
    "edge.length" => data.branch_lengths)

## cd /home/bkopper/projects/BDS_deterministic_map

res = calculate_tree_rates(data, model, Ds, Fs, Ss; verbose=false);
average_node_rates = res["average_node_rates"]

node = data.edges[:, 2]
lambda_avg = average_node_rates["λ"]
@rput nshifts
@rput phy
@rput node
@rput lambda_avg
@rput node_nshifts

R"""
library(ape)
library(ggtree)
library(tidytree)
library(ggplot2)
library(dplyr)
library(patchwork)
"""

R"""
class(phy) <- "phylo"
th <- max(node.depth.edgelength(phy))

df1 <- tibble("node" = 1:max(phy$edge),
            "nshifts" = node_nshifts)
#df1 <- rbind(df1, 
#            tibble("node" = length(phy$tip.label)+1,
#                   "nshifts" = 0))
x <- as_tibble(phy)

phydf <- merge(x, df1, by = "node")
df2 <- tibble("node" = 1:max(phy$edge),
                "lambda" = lambda_avg)
phydf <- merge(phydf, df2, by = "node")
td_phy <- as.treedata(phydf)

p1 <- ggtree(td_phy, aes(color = `nshifts`), size = 2) +
    geom_tiplab(size = 5) +
    theme(legend.position = c(0.2, 0.8)) +
    xlim(c(0.0, th + 10)) +
    ggtitle("Number of shifts") +
    scale_colour_gradient(low = "black", high = "red")

p2 <- ggtree(td_phy, aes(color = `lambda`),size =2) +
    geom_tiplab(size = 5) +
    theme(legend.position = c(0.2, 0.8)) +
    xlim(c(0.0, th + 10)) +
    ggtitle("Mean speciation rate")
p <- p1 + p2

ggsave("figures/cetacea_2state_nshifts.pdf", p, width = 300, height = 300, units = "mm")
""";



