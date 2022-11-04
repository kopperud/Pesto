using Diversification
using Distributions
using ProgressMeter
using CSV
using DataFrames
using StatsPlots

#### load simulated trees
fnames = Base.Filesystem.readdir("data/simulated_trees")

treefiles = ["data/simulated_trees/" * fname for fname in fnames]

data = Dict()
for treefile in treefiles
    phy = readtree(treefile)
    ntaxa = replace(split(treefile, "/")[3], ".tre" => "")
    ntaxa = parse(Int64, ntaxa)

    ρ = 1
    d = make_SSEdata(phy, "", ρ; include_traits = false)
    data[ntaxa] = d
end

## Set up the model
H = 0.587405
d1 = LogNormal(log(0.29), 3 * H)
d2 = LogNormal(log(0.20), 3 * H)

n = 6
speciation = make_quantiles(d1, n)
extinction = make_quantiles(d2, n)

k = n ^2
λ = zeros(k)
μ = zeros(k)
η = 0.1

for (i, (sp, ex)) in enumerate(Iterators.product(speciation, extinction))
    λ[i] = sp
    μ[i] = ex
end

model = SSEconstant(λ, μ, η)

## Compute the results
results = Dict()
times = Dict()
@showprogress for (ntaxa, d) in data
    t1 = time()
    Ds, Fs = backwards_forwards_pass(model, d; verbose = false)
    res = calculate_tree_rates(d, model, Ds, Fs; verbose = false);
    t2 = time()
    times[ntaxa] = t2 - t1

    results[ntaxa] = res
end

rdf = CSV.read("output/runtimes_r.csv", DataFrame)

jdf = DataFrame(time = collect(values(times)), ntaxa = collect(keys(times)))
sort!(jdf, :ntaxa)

ps = []
for i in 1:2
    p = plot(xlab = "ntaxa", ylab = "time (seconds)", legend = :bottomright, title = "Arithmetic")

    if i == 1
        plot!(p, xscale = :log10, yscale = :log10, title = "Log-scale")

        for i in 1:length(rdf[:,:ntaxa])
            annotate!(p, rdf[i, :ntaxa], rdf[i, :times], text(rdf[:,:ntaxa][i], 6, :left, :top))
        end
        for i in 1:length(jdf[:,:ntaxa])
            annotate!(p, jdf[i, :ntaxa], jdf[i, :time], text(jdf[:,:ntaxa][i], 6, :left, :top))
        end
    end
    plot!(p, jdf[:, :ntaxa], jdf[:, :time], label = "", color = "black")
    scatter!(p, jdf[:, :ntaxa], jdf[:, :time], label = "Julia", color = "black")

    ## R
    plot!(p, rdf[:, :ntaxa], rdf[:, :times], label = "", color = "orange")
    scatter!(p, rdf[:, :ntaxa], rdf[:, :times], label = "R", color = "orange")

    append!(ps, [p])
end

runtimeplot = plot(ps...)

savefig(runtimeplot, "figures/runtime_36states_100knotsperbranch.pdf")








