## q(t)

using Diversification
using Distributions
using QuadGK
using Plots
using ProgressMeter
using BenchmarkTools
using ForwardDiff
using RCall
using DataFrames
using CSV

phy = readtree(Diversification.path("primates.tre"))
ρ = 0.67
data = make_SSEdata2(phy, ρ)
λml, μml = estimate_constant_bdp(data)

## set up States
H = 0.587
n = 6

dλ = LogNormal(log(λml), H)
dμ = LogNormal(log(µml), H)

λ, μ = allpairwise(make_quantiles(dλ, n), make_quantiles(dμ, n))
η = Diversification.optimize_eta(λ, μ, data; lower = 0.00001, upper = 0.1)

model = SSEconstant(λ, μ, η[1])
res = bds(model, data)

function fq(model, F, E, K)
    v = (-1) .* ones(K)
    v[1] = 1

    w = (2-K) .* ones(K)
    w[1] = -1

    dF(t) = ForwardDiff.derivative(F, t)

    q(t) = ((v' * (dF(t) .- (2 .* model.λ .* E(t) .- model.λ .- model.μ) .* F(t))) / (2 .* w' * F(t)))[1]
    return(q)
end

E = Diversification.extinction_probability(model, data)

idx = 1

#q = fq(model, Fs[idx], E, 36)
#x = range(Fs[idx].t[end], Fs[idx].t[1], length = 100)
#plot(x, q.(x), linetype = [:line])

Ns = zeros(length(Fs))
@showprogress for (edge_idx, F) in Fs
    println(edge_idx)
    K = 36
    q = fq(model, F, E, K)
    a = F.t[end]
    b = F.t[1]

    N = quadgk(q, a, b, rtol = 1e-06)[1]
    Ns[edge_idx] = N
end
Ns_per_time = Ns ./ data.branch_lengths
h = plot(
    histogram(Ns, xlab = "N", title = "N per branch", lab = "branch", ylab = "frequency"), 
    histogram(Ns_per_time, xlab = "N/blength", title = "N per branch", lab = "branch", ylab = "frequency")
)
savefig(h, "figures/histogram_N_branch.pdf")


scplot = plot(data.branch_lengths, Ns, linetype = :scatter, 
xlab = "branch length", ylab = "E[N]", lab = "branch index")
savefig(scplot, "figures/scatter_N_branchlength")


## Reorder for ape indices
ancestors = Diversification.make_ancestors(data)
m = zeros(maximum(data.edges))
m2 = zeros(maximum(data.edges))

for i in 1:maximum(data.edges)
    if i == length(data.tiplab)+1
        m[i] = NaN
        m2[i] = NaN
    else
        edge_idx = ancestors[i]
        node_val = Ns[edge_idx]
        m[i] = node_val
        m2[i] = Ns_per_time[edge_idx]
    end
end



phy = res.phy

RCall.@rput phy
RCall.@rput m
RCall.@rput m2

RCall.R"""
class(phy) <- "phylo"
th <- max(ape::node.depth.edgelength(phy))

df1 <- tibble::tibble("node" = 1:max(phy$edge),
        "N" = m,
        "logN" = log(m),
        "N per time" = m2)
x <- tidytree::as_tibble(phy)

phydf <- merge(x, df1, by = "node")
td_phy <- tidytree::as.treedata(phydf)

p1a <- ggtree::ggtree(td_phy, ggplot2::aes(color = `N`))
    #ggplot2::scale_color_gradient(limits = c(0,5))
    #ggplot2::geom_text(ggplot2::aes(label=node))

ggplot2::ggsave("figures/number_of_transitions.pdf", p1a)
"""

## magnitude of change
res = bds(model, data)

Ds, Fs = backwards_forwards_pass(model, data);
Ps = ancestral_state_probabilities(data, model, Ds, Fs)

function mag(F, P, model)
    a = F.t[1]
    b = F.t[end]
    magnitude = sum(model.λ .* (P(b) - P(a)))
    return(magnitude)
end

mags = zeros(maximum(data.edges)-1)
for (edge_idx, p) in Fs
    m1 = mag(Fs[edge_idx], Ps[edge_idx], model)
    mags[edge_idx] = m1
end

deltasp_N = scatter(Ns, mags, xlab = "E[N]", ylab = "E[Δλ]", lab = "edge index")
savefig(deltasp_N, "figures/Nchanges_vs_delta_lambda.pdf")



## revbayes to ape indices
R"""
library(ape)
phy <- read.nexus("data/bears.tre")
source("scripts/matchnodes.R")
mn <- matchNodes(phy)
"""
@rget mn
mn[!,:Rev] = convert(Vector{Int64}, mn[!,:Rev])
Rev_to_R = Dict()
R_to_Rev = Dict()
for row in eachrow(mn)
    Rev_to_R[row["R"]] = row["Rev"]
    R_to_Rev[row["Rev"]] = row["R"]
end

## load the RevBayes log files
df = DataFrame(CSV.File("output/bears_scm.log"))
DataFrames.names(df)

using StatsBase

Nscm = zeros(14)
VarNscm = zeros(14)
for edge_idx in 1:14
    node_index = data.edges[edge_idx, 2]

    if node_index != 9
        Rev_index = R_to_Rev[node_index]
        Nscm[edge_idx] = mean(df[!,"num_shifts["*string(Rev_index)* "]"])
        VarNscm[edge_idx] = var(df[!,"num_shifts["*string(Rev_index)* "]"])
    end
end


p1 = plot(Nscm, Ns, linetype = :scatter, xlab = "Integrate q(t) over (t0, t1)", ylab = "RevBayes SCM")
p2 = plot(Nscm, Ns, yerror = VarNscm, linetype = :scatter)
plot(p1, p2)


