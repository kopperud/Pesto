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
using LinearAlgebra

phy = readtree("data/cetaceans.tre")
ρ = 1.0
data = make_SSEdata2(phy, ρ)

#λ = [0.1, 0.3, 0.1, 0.3]
#μ = [0.05, 0.10, 0.15, 0.20]
λ = [0.2, 0.1]
μ = [0.05, 0.1]
η = 0.1
K = length(λ)

model = SSEconstant(λ, μ, η)

res = birth_death_shift(model, data)

Ds, Fs = Diversification.backwards_forwards_pass(model, data);
Ps = Diversification.ancestral_state_probabilities(data, model, Ds, Fs)
E = extinction_probability(model, data)

## revbayes to ape indices
R"""
library(ape)
phy <- read.tree("data/cetaceans.tre")
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

## load the log file
## load the RevBayes log files
df = DataFrame(CSV.File("output/BDS_SCM_cetaceans_2.log"))

nbranches = length(data.branch_lengths)
Nscm = zeros(nbranches)
root_index = length(data.tiplab)+1
for edge_idx in 1:nbranches
    node_index = data.edges[edge_idx, 2]

    if node_index != root_index
        Rev_index = R_to_Rev[node_index]
        Nscm[edge_idx] = mean(df[!,"num_shifts["*string(Rev_index)* "]"])
    end
end


function my_ode1!(dN, N, p, t)
    F, D, K, dF, dG, η = p

    Ft = F(t)
    sumF = sum(Ft)
    Dt = D(t)
    sumD = sum(Dt)

    for i in 1:K
        dN[i] = 0.0
        for j in 1:K
            if j != i
                dN[i] -= (Ft[i] / sumF) * (η/(K-1)) * (Dt[j] / sumD)
            end
        end
    end
end

function my_ode2!(dN, N, p, t)
    F, D, K, dF, dG, η = p

    dN[:] = - (η/(K-1)) .* (1.0/sum(F(t))) .* (1.0 / sum(D(t))) .* 
                            (1 .- LinearAlgebra.I(K)) .* (F(t) * D(t)') * ones(K)
end

function my_ode3!(dN, N, p, t)
    F, D, K, dF, dG, η = p

    dN[:] = - (η/(K-1)) .* (1.0 / sum(D(t))) .* (1 .- LinearAlgebra.I(K)) * D(t)
end

function my_ode4!(dN, N, p, t)
    F, D, K, dF, dG, η = p

    dN[:] = ones(K) .* -η ./ K
end

function my_ode5!(dN, N, p, t)
    F, D, K, dF, dG, η = p

    dN[:] = - dF(t) ./ sum(dF(t))
end

function my_ode6!(dN, N, p, t)
    F, D, K, dF, dG, η = p

    dN[:] = - (η/(K-1)) .* sum((1 .- I(K)) .* ((dF(t) ./ (D(t) .* sum(dF(t)))) * D(t)'), dims = 2) 
end

function my_ode7!(dN, N, p, t)
    F, D, K, dF, dG, η = p

    trprob = (dF(t) * F(t)') ./ sum(dF(t) * F(t)', dims = 1)

    dN[:] = - (η/(K-1)) .* sum((1 .- I(K)) .* trprob, dims = 2)
end

function my_ode8!(dN, N, p, t)
    F, D, K, dF, dG, η = p

    #denom = sum(dG(t))
    denom = 1.0

    dN[:] = -(η/(K-1)) .* sum((1 .* I(K)) .* (dG(t) * ones(K)') ./ denom, dims = 2)
end

ForwardDiff.derivative(Fs[1], 12)
ForwardDiff.derivative(Fs[1], 12) ./ sum(ForwardDiff.derivative(Fs[1], 12))

my_odes = [
    my_ode1!,
    my_ode2!,
    my_ode3!,
    my_ode4!,
    my_ode5!,
    my_ode6!,
    my_ode7!,
    my_ode8!
]

S = zeros(length(Fs), length(my_odes))
for (i, my_ode) in enumerate(my_odes)
    my_ode = my_odes[i]

    for edge_idx in 1:nbranches
        a = Fs[edge_idx].t[1]
        b = Fs[edge_idx].t[end]

        tspan = (a, b)
        derivF(t) = ForwardDiff.derivative(Fs[edge_idx], t)
        G(t) = Fs[edge_idx](t) ./ sum(Fs[edge_idx](t))
        dG(t) = ForwardDiff.derivative(G, t) 


        p = [Fs[edge_idx], Ds[edge_idx], K, derivF, dG, η]
        N0 = zeros(K)
        prob = ODEProblem(my_ode, N0, tspan, p)

        sol = solve(prob)
        S[edge_idx, i] = sum(sol[end])
    end
end

ps5 = []
for i in 1:length(my_odes)
    p = plot(S[:,i], Nscm, linetype = :scatter, xlab = "ODE", ylab = "SCM", label = "#state changes")
    plot!(p, [0.0, 0.2], [0.0, 0.2], linestyle = :dash, label = "One-to-one")
    append!(ps5, [p])
end
plot(ps5..., size = (700, 700))


function Amatrix(t)
    Q = -I(K) .* η .+ (1 .- I(K)) .* (η/(K-1))
    res = (λ .+ μ .- 2 .* λ .* E(t)) .* I(K) .- Q
    return(res)
end

Amatrix(30) * Fs[1](30)

ForwardDiff.derivative(Fs[1], 30)

Δt = -0.01
first = Amatrix(30) * [0.0, 1.0] .* Δt .+ [0.0, 1.0]
second = Amatrix(30) * [1.0, 0.0] .* Δt .+ [1.0, 0.0]


function trans_prob!(dP, P, p, t)
    K, E, η = p

    Q = -I(K) .* η .+ (1 .- I(K)) .* (η/(K-1))
    A = diagm(λ .+ μ .- 2 .* λ .* E(t)) .- Q

    #denom = sum(A * I(K), dims = 1)
    dP[:] = - A * I(K) * ones(K)
end

#P0 = (1/K) .* ones(K) * ones(K)'
P0 = [
    0.0, 0.0
]

p = [K, E, η]
edge_idx = 104
tspan = (Fs[edge_idx].t[1], Fs[edge_idx].t[end])

prob = ODEProblem(trans_prob!, P0, tspan, p)
sol = solve(prob)

plot(sol, linetype = [:line])

U, S, V = svd(Amatrix(30))

svd([
    -1.0 1/2 1/2
    1/2 -1.0 1/2
    1/2 1/2 -1.0
])

U * diagm(S) * V

function dNbinary(dN, N, p, t)
    a, b, G, D = p

    dN[1] = a(t) * G(t)[1]#* (D(t)[2] ./ sum(D(t)))
    dN[2] = b(t) * G(t)[2]#* (D(t)[1] ./ sum(D(t)))
end

n_analytical2 = zeros(nbranches)
for edge_idx in 1:nbranches
    G(t) = Fs[edge_idx](t) ./ sum(Fs[edge_idx](t))
    dG(t) = ForwardDiff.derivative(G, t)

    #af(t) = ForwardDiff.derivative(Fs[edge_idx], t)[1] / (Fs[edge_idx](t)[2] - Fs[edge_idx](t)[1])
    #bf(t) = ForwardDiff.derivative(Fs[edge_idx], t)[2] / (Fs[edge_idx](t)[1] - Fs[edge_idx](t)[2])
    af(t) = ForwardDiff.derivative(G, t)[1] / (G(t)[2] - G(t)[1])
    bf(t) = ForwardDiff.derivative(G, t)[2] / (G(t)[1] - G(t)[2])

    p = [af, bf, G, Ds[edge_idx]]
    tspan = (Fs[edge_idx].t[1], Fs[edge_idx].t[end])

    N0 = zeros(K)
    prob = ODEProblem(dNbinary, N0, tspan, p)

    sol = solve(prob)
    n_analytical2[edge_idx] = sum(sol[end])
end

edge_idx = 104
af(t) = ForwardDiff.derivative(Fs[edge_idx], t)[1] / (Fs[edge_idx](t)[2] - Fs[edge_idx](t)[1])
bf(t) = ForwardDiff.derivative(Fs[edge_idx], t)[2] / (Fs[edge_idx](t)[1] - Fs[edge_idx](t)[2])
x = collect(range(Fs[edge_idx].t[1], Fs[edge_idx].t[end], length = 50))
plot(x, af.(x))
plot!(x, bf.(x))

plot(n_analytical2, Nscm, linetype = :scatter, xlab = "analytical", ylab = "SCM")
















ForwardDiff.derivative(Fs[1], 30) ./ (M * Fs[1](30))












Nscm[1]





















hcat(S, η .* data.branch_lengths)


## posterior transition probability

function dN10(dP, P, p, t)
    F, D, K, dF, dG, η = p

    denom = sum(dG(t))

    dP[:] = -(η/(K-1)) .* sum((1 .* I(K)) .* (dG(t) * ones(K)') ./ denom, dims = 2)
end

edge_idx = 1
a = Fs[edge_idx].t[1]
b = Fs[edge_idx].t[end]

tspan = (a, b)
derivF(t) = ForwardDiff.derivative(Fs[edge_idx], t)
trprob(t) = (derivF(t) * D(t)') ./ sum(derivF(t) * D(t)')
G(t) = Fs[edge_idx](t) ./ sum(Fs[edge_idx])
Gderiv(t) = ForwardDiff.derivative(G, t)

p = [Fs[edge_idx], Ds[edge_idx], K, derivF, Gderiv, η]
N0 = zeros(K)
prob = ODEProblem(dN10, N0, tspan, p)
sol = solve(prob)

plot(sol)

sol[end] |> sum

Nscm
Nscm[1]

P = derivF(a) * Ds[1](a)'
P = P ./ sum(P, dims = 1)


## load cetaceans

df1 = DataFrame(CSV.File("output/bears_scm_run_1.log"))
df2 = DataFrame(CSV.File("output/bears_scm_run_2.log"))
df = vcat(df1, df2)



















