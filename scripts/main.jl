using CSV
using RCall
using DifferentialEquations
using Distributions

include("datatypes.jl")
include("utils.jl")
include("diversitree_bisse.jl")
include("postorder.jl")
include("preorder.jl")
include("ODE.jl")
include("logLroot.jl")

function make_quantiles(d, k)
    n_breakpoints = k + 2
    breakpoints = range(0.0, 1.0, length = n_breakpoints)
    quantiles = quantile.(d, breakpoints)[2:end-1]
end

if false
    λ = [2.0, 1.0]
    μ = [0.5, 0.1]
else
#    λ = [0.1, 0.5, 1.0, 1.5, 2.0]
#    μ = [0.01, 0.2, 0.5, 1.0, 1.5]
    d1 = LogNormal(log(0.5))
    d2 = LogNormal(log(0.3))
    k = 10

    λ = make_quantiles(d1, k)
    μ = make_quantiles(d2, k)
end
η = 0.05
ρ = 1.0
k = length(λ)

model = SSEconstant(λ, μ, η)

treefile = "data/bears.tre"
datafile = "data/bears.csv"

if false
    println("Diversitree-BISSE model:\n")
    asr = anc_state_prob_bisse(treefile, datafile, model)
    println("Ancestral state marginal probabilities:")
    display(asr)
end

println("\t\t")
println("Backwards-forwards method:")
phy = readtree(treefile)
data = make_SSEdata(phy, datafile, ρ; include_traits = false)
D_ends, Ds, sf, E = postorder(model, data)
logL = logL_root(model, data)
println("logL: \t", logL)

ASP, Fs = preorder(model, data, E, D_ends)
display(ASP)
println()

Ps = Dict()
for edge_idx in 1:(maximum(data.edges)-1)
    ## Try to normalize
#    Dn(t) = D(t) ./ sum(D(t))
#    Fn(t) = F(t) ./ sum(F(t))
#    P(t) = Fn(t) .* Dn(t) ./ (sum(Fn(t) .* Dn(t)))
    Ps[edge_idx] = t -> Fs[edge_idx](t) .* Ds[edge_idx](t) ./ (sum(Fs[edge_idx](t) .* Ds[edge_idx](t)))
end

average_branch_rates = Dict()
for (key, P) in Ps
    times = Fs[key].t
    average_branch_rates[key] = mean([sum(λ .* P) for P in Ps[key].(times)])
end

average_node_rates = zeros(maximum(data.edges))
for i in 1:maximum(data.edges)
    if i == length(data.tiplab)+1
        average_node_rates[i] = NaN      
    else
        edge_idx = findall(data.edges[:,2] .== i)[1]
        node_val = average_branch_rates[edge_idx]
        average_node_rates[i] = node_val      
    end
end

@rput average_node_rates
phy = Dict("edge" => data.edges,
          "tip.label" => data.tiplab,
          "Nnode" => length(data.tiplab)-1,
         "edge.length" => data.branch_lengths)
@rput phy
R"""
library(ape)
library(ggtree)
library(tidytree)

class(phy) <- "phylo"

d <- tibble("node" = 1:max(phy$edge),
            "Speciation rate" = average_node_rates)
x <- as_tibble(phy)

phydf <- merge(x, d, by = "node")
phy <- as.treedata(phydf)

p <- ggtree(phy, aes(color = `Speciation rate`)) +
  geom_tiplab()
"""

using StatsPlots

i = 13
#p1 = plot(Fs[i].t, [mean(Ps[i](t) .* model.λ) for t in Fs[i].t], xflip = true, xlab = "time (Ma)", ylab = "P(t) for edge "*string(i))
#p2 = plot(times1, [mean(Ps[i](t) .* model.λ) for t in times1], xflip = true, xlab = "time (Ma)", ylab = "P(t) for edge "*string(i))

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


