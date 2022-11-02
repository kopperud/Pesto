using RCall
using Distributions
using StatsPlots
using Diversification
using DataFrames
using CSV
using ProgressMeter

##############################
##
##   Load the data files
##
###############################
treefiles = []

filenames = Base.Filesystem.readdir("data/condamine_2019")
treefiles = filenames[contains.(filenames, r".tre")]
treefiles = [fn for fn in treefiles if !contains(fn, r"Picidae2")]
samplingfiles = filenames[contains.(filenames, r".Rev")]

data = []
for (treefile, samplingfile) in zip(treefiles, samplingfiles)
    phy = readtree("data/condamine_2019/" * treefile)
    
    txt = readlines("data/condamine_2019/" * samplingfile)[1]
    ρ = parse(Float64, split(txt, " ")[4])
 
    x = make_SSEdata(phy, "", ρ; include_traits = false)
    append!(data, [x])
end

##############################
##
##   Set up the model, common for all datasets
##
###############################
H = 0.587405
d1 = LogNormal(log(0.08251888), 3 * H)
d2 = LogNormal(log(0.08251888/2.0), 3 * H)

n = 6
speciation = make_quantiles(d1, n)
extinction = make_quantiles(d2, n)
η = 0.1

k = n ^2
λ = zeros(k)
μ = zeros(k)

for (i, (sp, ex)) in enumerate(Iterators.product(speciation, extinction))
    λ[i] = sp
    μ[i] = ex
end

model = SSEconstant(λ, μ, η)

## Calculate the backwards-forwards pass equations
phys = Dict()
avg_branch_rates = Dict()
cladenames = [replace(split(x, "_")[4], ".tre" => "") for x in treefiles]

@showprogress for (d, cladename) in zip(data, cladenames)
    Ds, Fs = backwards_forwards_pass(model, d; verbose = false) 
    res = calculate_tree_rates(d, model, Ds, Fs; verbose = false);
    average_edge_rates = res["average_node_rates"]
    avg_branch_rates[cladename] = average_edge_rates

    phy = Dict("edge" => d.edges,
          "tip.label" => d.tiplab,
          "Nnode" => length(d.tiplab)-1,
         "edge.length" => d.branch_lengths)
    phys[cladename] = phy
end

hs = Dict()
for (i, cladename) in enumerate(cladenames)
    br = avg_branch_rates[cladename]

    if i == 1
        labels = ["λ", "μ"]
    else
        labels = ["", ""]
    end
    maxlambda = maximum(br["λ"][.!isnan.(br["λ"])])
    maxmu = maximum(br["μ"][.!isnan.(br["μ"])])
    h = histogram(br["λ"], label = labels[1], bins = range(0.0, maxlambda, length = 15), title = cladename, xlab = "Rate", xrotation = 90, xlim = (0.0, 0.6))
    histogram!(h, br["μ"], bins = range(0.0, maxmu, length = 10), label = labels[2])

    #, xlim = (0.0, 1.0)
    hs[cladename] = h    
end

big_histogram = plot(collect(values(hs))[1:20]..., size = (800, 800))

savefig(big_histogram, "figures/big_histogram.pdf")

f(λ, μ) = logpdf(d1, λ) + logpdf(d2, μ)
surface(0.001:0.01:maximum(λ), 0.001:0.01:maximum(μ), f, xlab = "λ", ylab = "μ", zlab = "logL")


@showprogress for cladename in cladenames
    lambda_average = avg_branch_rates[cladename]["λ"]
    phy = phys[cladename]

    @rput lambda_average
    @rput phy
    @rput cladename
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

    fname <- paste0("figures/condamine_2019/", cladename, ".pdf")
    ggsave(fname, p1a, width = 160, height = 160, units = "mm")
    """
end

