using Diversification
using Distributions
using Plots
using ProgressMeter
using BenchmarkTools
using ForwardDiff
using RCall
using DataFrames
using CSV
using StatsBase
using DifferentialEquations
using LaTeXStrings
using Measures

phy = readtree(Diversification.path("bears.tre"))
ρ = 1.0
data = make_SSEdata2(phy, ρ)

λ = [0.1, 0.3, 0.1, 0.3]
μ = [0.05, 0.10, 0.15, 0.20]
η = 0.1

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
    Rev_to_R[row["Rev"]] = row["R"]
    R_to_Rev[row["R"]] = row["Rev"]
end

## load the RevBayes log files
df1 = DataFrame(CSV.File("output/bears_scm_run_1.log"))
df2 = DataFrame(CSV.File("output/bears_scm_run_2.log"))
df = vcat(df1, df2)
#df = df[1:20000, :]
DataFrames.names(df)

Nscm = zeros(14)
for edge_idx in 1:14
    node_index = data.edges[edge_idx, 2]

    if node_index != 9
        Rev_index = R_to_Rev[node_index]
        Nscm[edge_idx] = mean(df[!,"num_shifts["*string(Rev_index)* "]"])
    end
end

S = data.branch_lengths .* η

std_err_mean = sqrt.(Nscm/size(df)[1])

is_tip = data.edges[:,2] .< 9
scm_s = plot(S[is_tip], Nscm[is_tip], linetype = :scatter,
            ylab = "RevBayes\nStochastic Character Maps", yerror = std_err_mean, label = "S (tip edge)",
            xlab = L"Present approach, $S(t) = (t_a - t)η$", color = "black", grid = :false)
plot!(scm_s, S[.!is_tip], Nscm[.!is_tip], linetype = :scatter, color = "orange", label = "S (internal edge)")
linemax = maximum([maximum(Nscm), maximum(S)])
plot!(scm_s, [0.0, linemax], [0.0, linemax], label = "One-to-one", color = "black")        

savefig(scm_s, "figures/ode_vs_scm.pdf")

## How many iterations do I need?
#df3 = DataFrame(CSV.File("output/bears_scm.log"))
ns = [1000, 2000, 5000, 10000, 15000, 20000, 50000, 100000]

ys = zeros(length(ns), 14)
for (i,n) in enumerate(ns)
    df3 = df[1:n,:]

    for j in 1:14
        node_index = data.edges[j, 2]
        Rev_index = R_to_Rev[node_index]

        y = mean(df3[!,"num_shifts["*string(Rev_index)* "]"])
        ys[i,j] = y
    end
end

iterplot = plot(xlab = "samples (iid)", ylab = "mean(N)")

## time at tip
for edge_idx in 1:14
    if is_tip[edge_idx]
        tb = Fs[edge_idx].t[end]
        println("edge $edge_idx: t = $tb")
    end
end


iterplots = []
for edge_index in 1:14
    myplot = plot(ns, ys[:,edge_index], linetype = [:scatter, :line], color = "black", label = "")#, ylim = (0.0, 0.005))
    append!(iterplots, [myplot])
end
ps2 = plot(iterplots..., xrotation = 90, size = (600, 600))

savefig(ps2, "figures/how_many_iterations_scm.pdf")


### Try to simulate under a branch
edge_idx = 1

## In this time span
a = Fs[edge_idx].t[1]
b = Fs[edge_idx].t[end]

# How many time steps?
ntimes = 20000
Δt = (b-a) / (ntimes -1)
times = collect(range(a, b, length = ntimes))


# shift probability in Δt
η = 0.1
shift_probability = -η * Δt

# initial value probability
starting_distribution = Categorical(Ps[edge_idx](a))

y = zeros(Int64, ntimes, K)
nshifts = 0
n_iters = 100000

@showprogress for i in 1:n_iters
    state = rand(starting_distribution)
    y[1, state] += 1
    for j in 2:ntimes
        u = rand(1)[1]
        if u < shift_probability
            p = Ps[edge_idx](times[j]) ## should this be j-1 ?
            ## rescale so that current state prob is 0
            p[state] = 0.0
            p = p ./ sum(p)
            ## redraw from the categorical Distributions
            current_distribution = Categorical(p)

            state = rand(current_distribution)[1]
            nshifts += 1
        end
        y[j,state] += 1
    end
end
mean(y, dims = 1)
nshifts / n_iters

## 1.07643 with P(t[j-1]) and 100k iters
## 1.07923 with P(t[j]) and 100k iters

Ns_ode[4,1]

Ps[edge_idx](a)
nshifts



























