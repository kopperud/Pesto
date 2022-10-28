using DataFrames, CSV, StatsPlots
using Diversification
using RCall

## Read RevBayes analytical expression (Euler's method)

filenames = Base.Filesystem.readdir("output/ASP_Euler")
#nts = range(100, 10100; step = 1000)

ds = Dict("bears" => [],
                "primates" => [])

for filename in filenames
    s = split(filename, "_")
    dname = s[3]
    nTimeSlices = parse(Int64, split(s[4], ".")[1])
#    push!(nts, nTimeSlices)

    df = CSV.read("output/ASP_Euler/" * filename, DataFrame)
    df[:, :nTimeSlices] = [nTimeSlices]
    ## sort by nTimeSlices
    #df = sort(df, :nTimeSlices)
    append!(ds[dname], [df])
end

datasets = Dict(key => vcat(val...) for (key, val) in ds)



#df = datasets["bears"]
#dataset = "bears"

for (dataset, df) in datasets
    sort!(df, :nTimeSlices)
    ## DO THE JULIA IMPLEMENTATION OF BDS
    treefile = "data/" * dataset * ".tre"
    phy = readtree(treefile)
    ρ = 1.0

    data = make_SSEdata(phy, "", ρ; include_traits = false)
    H = 0.587405
    ntaxa = length(data.tiplab)
    tree_length = sum(data.branch_lengths)
    rate_mean = (ntaxa - 2) / tree_length 
    d1 = LogNormal(log(rate_mean), 3 * H)
    d2 = LogNormal(log(rate_mean/2.0), 3 * H)

    k = 20
    λ = make_quantiles(d1, k)
    μ = make_quantiles(d2, k)
    η = 0.1

    model = SSEconstant(λ, μ, η)

    Ds, Fs = backwards_forwards_pass(model, data; verbose = false)
    res = calculate_tree_rates(data, model, Ds, Fs; verbose = false);
    average_node_rates = res["average_node_rates"]

    phy = Dict("edge" => data.edges,
              "tip.label" => data.tiplab,
              "Nnode" => length(data.tiplab)-1,
              "edge.length" => data.branch_lengths)

    R"""
    source("scripts/matchnodes.R")    

    mn <- matchNodes(phy)
    """
    @rget mn
    mn = convert.(Int64, mn)

    y = Matrix(df[!, names(df, r"avg_lambda")])

    plots = []
    nTimeSlices = df[!, "nTimeSlices"]

    for (i, RevBayes_index) in enumerate(mn[!, "Rev"])
        p = plot(nTimeSlices, y[:,RevBayes_index], title = "branch "* string(i), lab = "RB Analytical", xrotation = 90, legend = :bottomright)

        if i > 1
            plot!(p, legend = :none)
        end
        if i == 1
            plot!(p, ylab = "mean branch λ", xlab = "nTimeSlices")
        end

        ## add JL rates
        plot!(p, [nTimeSlices[1], nTimeSlices[end]], [average_node_rates["λ"][i], average_node_rates["λ"][i]], lab = "jl impl")


        append!(plots, [p])
    end



    if dataset == "bears"
        myplot = plot(plots[1:14]..., size = (800, 800))
    elseif dataset == "primates"
        for i in 1:length(plots)
            if i > 10
                plot!(plots[i], xticks = false)
            end
        end
        myplot = plot(plots[1:400]..., size = (1500, 8000), layout = (40, 10))
#        myplot = plot(plots[1:100]..., size = (2000, 2000))
    end

    savefig(myplot, "figures/euler_rb_" * dataset * ".pdf")
end




