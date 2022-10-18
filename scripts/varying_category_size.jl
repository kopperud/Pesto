using Diversification
using Distributions
using BenchmarkTools
using StatsPlots

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

bs = [] ## benchmarks
ks = [2, 4, 8, 16, 32, 64, 128, 256]
for k in ks 
    λ = make_quantiles(d1, k)
    μ = make_quantiles(d2, k)
    η = 0.1
    k = length(λ)

    model = SSEconstant(λ, μ, η)

    ## Calculate the backwards-forwards pass equations
#    display(model)
    b = @benchmark backwards_forwards_pass($model, data; verbose = false) samples = 50
#    run(b, samples = 50)
    Ds, Fs = backwards_forwards_pass(model, data; verbose = false) 
    res = calculate_tree_rates(data, model, Ds, Fs; verbose = false);
    average_node_rates = res["average_node_rates"]

    append!(bs, [b])
    print(".")
end

## Plot the result of how number of categories scales with time
tplot = plot()
for (k, b) in zip(ks, bs)
    violin!(tplot, [k], [b.times .* 10^(-9)], lab = string(k)*" categories")
end
plot!(tplot, [ks[1], ks[end]], [mean(bs[1].times), mean(bs[end].times)] .* 10^(-9), lab = "Linear scaling")
plot!(tplot, legend = :topleft, ylab = "seconds (50 replicates)")

## Save the plot
savefig(tplot, "figures/varying_number_categories.pdf")
