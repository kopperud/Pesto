using CSV
using RCall
using DifferentialEquations
using Distributions
using ProgressMeter

include("datatypes.jl")
include("utils.jl")
include("diversitree_bisse.jl")
include("postorder.jl")
include("preorder.jl")
include("ODE.jl")
include("logLroot.jl")
include("constantBDP.jl")

#function make_quantiles(d, k)
#    n_breakpoints = k + 2
#    breakpoints = range(0.0, 1.0, length = n_breakpoints)
#    quantiles = quantile.(d, breakpoints)[2:end-1]
#end

function make_quantiles(d, k)
    quantiles = zeros(k)
    for i in 1:k
        p = (i-0.5)/k
        quantiles[i] = quantile(d, p)
    end
    return(quantiles)
end

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
if false
    λ = [2.0, 1.0]
    μ = [0.5, 0.1]
else
    H = 0.587405
    d1 = LogNormal(log(0.08251888), 3 * H)
    d2 = LogNormal(log(0.08251888/2.0), 3 * H)

    k = 20
    λ = make_quantiles(d1, k)
    μ = make_quantiles(d2, k)
end
η = 0.1
k = length(λ)

model = SSEconstant(λ, μ, η)

if false
    println("Diversitree-BISSE model:\n")
    asr = anc_state_prob_bisse(treefile, datafile, model)
    println("Ancestral state marginal probabilities:")
    display(asr)
end

println("\t\t")
println("Backwards-forwards method:")
D_ends, Ds, sf, E = postorder(model, data)
#logL = logL_root(model, data)
#println("logL: \t", logL)

ASP, Fs = preorder(model, data, E, D_ends)
display(ASP)
println()

println("Calculating state probabilities")
Ps = Dict()
@showprogress for edge_idx in 1:(maximum(data.edges)-1)
   Ps[edge_idx] = t -> Fs[edge_idx](t) .* Ds[edge_idx](t) ./ (sum(Fs[edge_idx](t) .* Ds[edge_idx](t)))
end

println("Calculating average branch rates")
average_branch_rates = Dict()
for (rate, rate_name) in zip([λ, μ], ("λ", "μ"))
    d = Dict()
    @showprogress for (key, P) in Ps
        times = Fs[key].t
        d[key] = mean([sum(rate .* P) for P in Ps[key].(times)])
    end
    average_branch_rates[rate_name] = d
end

println("Reordering for ape node indices")
average_node_rates = Dict()
for (rate, rate_name) in zip([λ, μ], ("λ", "μ"))
    m = zeros(maximum(data.edges))
    @showprogress for i in 1:maximum(data.edges)
        if i == length(data.tiplab)+1
            m[i] = NaN      
        else
            edge_idx = findall(data.edges[:,2] .== i)[1]
            node_val = average_branch_rates[rate_name][edge_idx]
            m[i] = node_val      
        end
    end
    average_node_rates[rate_name] = m
end

lambda_average = average_node_rates["λ"]
@rput lambda_average
phy = Dict("edge" => data.edges,
          "tip.label" => data.tiplab,
          "Nnode" => length(data.tiplab)-1,
         "edge.length" => data.branch_lengths)
@rput phy
R"""
library(ape)
library(ggtree)
library(tidytree)
library(ggplot2)

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
rates <- select(df, starts_with("avg_lambda"))

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
plot!(vplot, ylim = (0.0, 0.5))


using StatsPlots

ps = []
for i in 1:14
    times1 = range(minimum(Fs[i].t), maximum(Fs[i].t); length = 100)
    p1 = plot(times1, hcat(Ds[i].(times1)...)', title = "D")
    p2 = plot(times1, hcat(Fs[i].(times1)...)', title = "F")
    p3 = plot(times1, hcat(Ps[i].(times1)...)', title = "P")
    p = plot(p1, p2, p3, layout = (1,3), xflip = true, nrow = 1 )
    append!(ps, [p])
end

p = plot(d1); for q in λ plot!(p, [q, q], [0.0, exp(loglikelihood(d1, q))], color = "green") end; plot(p)


