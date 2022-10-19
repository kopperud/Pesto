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
H = 0.587405

#tree_length = sum(data.branch_lengths)
#rate_mean = (num_total_species - 2.0) / tree_length
#rate_mean = rate_mean / 5.0
#d1 = LogNormal(log(rate_mean), 3 * H)
#d2 = LogNormal(log(rate_mean), 3 * H)
d1 = LogNormal(log(0.29), 3 * H)
d2 = LogNormal(log(0.20), 3 * H)

n = 6
speciation = make_quantiles(d1, n)
extinction = make_quantiles(d2, n)

k = n ^2
λ = zeros(k)
μ = zeros(k)

for (i, (sp, ex)) in enumerate(Iterators.product(speciation, extinction))
    λ[i] = sp
    μ[i] = ex
end

models = []
ηs = [0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01, 0.005, 0.001]
ηs = [0.04, 0.03, 0.02, 0.01, 0.005, 0.001]

for η in ηs 
    model = SSEconstant(λ, μ, η)
    append!(models, [model])
end

R"""
eta_plots <- list()
"""

l = []
for (i, model) in enumerate(models)
    ## Calculate the backwards-forwards pass equations
    Ds, Fs = backwards_forwards_pass(model, data; verbose = true) 
    res = calculate_tree_rates(data, model, Ds, Fs; verbose = false);

    average_node_rates = res["average_node_rates"]

    phy = Dict("edge" => data.edges,
          "tip.label" => data.tiplab,
          "Nnode" => length(data.tiplab)-1,
         "edge.length" => data.branch_lengths)

    lambda_average = average_node_rates["λ"]
    append!(l, [lambda_average])
    eta = model.η

    @rput lambda_average
    @rput phy
    @rput i
    @rput eta
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
        ggtitle(paste0("eta = ", eta)) +
        geom_tiplab(size = 1) +
        theme(legend.position = c(0.2, 0.8)) +
        xlim(c(0.0, th + 10)) +
    scale_color_gradient(limits = c(0.091, 0.5))

    eta_plots[[i]] <- p1a
    """;
end

println("min: ", minimum([minimum(x[x .> 0.0]) for x in l]))
println("max: ", maximum([maximum(x[x .> 0.0]) for x in l]))

R"""
p_combined <- Reduce("+", eta_plots)
ggsave("figures/primates_vaying_eta.pdf", p_combined, units = "mm", width = 400, height = 350)
"""
