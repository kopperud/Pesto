using RCall
using Distributions
using StatsPlots
using Diversification
using ProgressMeter

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
##   Set up the model variation
##
###############################

## Number of variation in hyperparameters

nhyper = 7 
model_sets = Dict("λ" => [], "μ" => [], "η" => [])

foo(x,i) = x * 2.0 ^(i-4)

## Central estimates
λhat = 0.29
μhat = 0.20
ηhat = 0.1
n = 30
H = 0.587405


## λ
for j in 1:nhyper
    λh = foo(λhat, j)
    η = ηhat

    d1 = LogNormal(log(λh), 3 * H)
    d2 = LogNormal(log(μhat), 3 * H)
    speciation = make_quantiles(d1, n)
    extinction = make_quantiles(d2, n)

    k = n^2
    λ = zeros(k)
    μ = zeros(k)

    for (i, (sp, ex)) in enumerate(Iterators.product(speciation, extinction))
    λ[i] = sp
    μ[i] = ex
    end

    model = SSEconstant(λ, μ, η)
    append!(model_sets["λ"], [model])
end

## μ
for j in 1:nhyper
    λh = λhat
    μh = foo(μhat, j)
    η = ηhat

    d1 = LogNormal(log(λhat), 3 * H)
    d2 = LogNormal(log(μh), 3 * H)
    speciation = make_quantiles(d1, n)
    extinction = make_quantiles(d2, n)

    k = n^2
    λ = zeros(k)
    μ = zeros(k)

    for (i, (sp, ex)) in enumerate(Iterators.product(speciation, extinction))
    λ[i] = sp
    μ[i] = ex
    end

    model = SSEconstant(λ, μ, η)
    append!(model_sets["μ"], [model])
end

## η
for j in 1:nhyper
    λh = λhat
    μh = μhat
    η = foo(ηhat, j)

    d1 = LogNormal(log(λh), 3 * H)
    d2 = LogNormal(log(μhat), 3 * H)
    speciation = make_quantiles(d1, n)
    extinction = make_quantiles(d2, n)

    k = n^2
    λ = zeros(k)
    μ = zeros(k)

    for (i, (sp, ex)) in enumerate(Iterators.product(speciation, extinction))
    λ[i] = sp
    μ[i] = ex
    end

    model = SSEconstant(λ, μ, η)
    append!(model_sets["η"], [model])
end

results = Dict()
for (rate_name, models) in model_sets
    l = []
    @showprogress for model in models
        ## Calculate the backwards-forwards pass equations
        Ds, Fs = backwards_forwards_pass(model, data; verbose = false) 
        res = calculate_tree_rates(data, model, Ds, Fs; verbose = false);
        append!(l, [res])
    end
    results[rate_name] = l
end

rmnan(x) = x[.!isnan.(x)]

ps = []

for rate_name in keys(model_sets)
    p = plot(xlim = (0.0, 8.0), legend = :topleft, xlabel = "Hyperparameter on " * rate_name, ylab = "Average branch rates")

    for (i, res) in enumerate(results[rate_name])
        λ = rmnan(res["average_node_rates"]["λ"])
        μ = rmnan(res["average_node_rates"]["μ"])
        if i == 1
            λlabel = "λ"
            μlabel = "μ"
        else
            λlabel = ""
            μlabel = ""
        end
        violin!(p, [i+0.05], λ, label = λlabel, color = "lightblue")
        violin!(p, [i-0.05], μ, label = μlabel, color = "red")
    end

    if rate_name == "λ"
        vals = λhat
    elseif rate_name == "μ"
        vals = μhat
    else
        vals = ηhat
    end

    xticks!([1:1:7;], [string(x) for x in foo.(vals, 1:7)])
    plot!(p, ylim = (0.0, 1.0))
    append!(ps, [p])
end

pcombined = plot(ps..., layout = (3, 1), size = (500, 700))

savefig(pcombined, "figures/hyperparameters.pdf")


#average_node_rates = res["average_node_rates"]
#
#phy = Dict("edge" => data.edges,
#      "tip.label" => data.tiplab,
#      "Nnode" => length(data.tiplab)-1,
#     "edge.length" => data.branch_lengths)
#
#lambda_average = average_node_rates["λ"]

#@rput lambda_average
#@rput phy
#R"""
#library(ape)
#library(ggtree)
#library(tidytree)
#library(ggplot2)
#library(dplyr)
#
#class(phy) <- "phylo"
#th <- max(node.depth.edgelength(phy))
#
#df1 <- tibble("node" = 1:max(phy$edge),
#            "Speciation rate" = lambda_average)
#x <- as_tibble(phy)
#
#phydf <- merge(x, df1, by = "node")
#td_phy <- as.treedata(phydf)
#
#p1a <- ggtree(td_phy, aes(color = `Speciation rate`)) +
#    geom_tiplab(size = 8) +
#    theme(legend.position = c(0.2, 0.8)) +
#    xlim(c(0.0, th + 10)) 
#""";

