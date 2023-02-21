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
using QuadGK

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
Ss = Diversification.ancestral_state_probabilities(data, model, Ds, Fs)
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

nbranches = length(data.branch_lengths)
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
    res = diagm(- (λ .+ μ .- 2 .* λ .* E(t))) .+ Q
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

### Attempt3

## ηᵢ(t) = \frac{2 λᵢ Eᵢ(t) Fᵢ(t) + dFᵢ(t)/dt - (λᵢ + μᵢ)Fᵢ(t)}{Fᵢ(t) - (1/(K-1))∑_{j ̸= i} Fⱼ(t)}

edge_idx = 104
Fs[edge_idx]
fderiv(t) = ForwardDiff.derivative(Fs[edge_idx], t)
pderiv(t) = ForwardDiff.derivative(Ps[edge_idx], t)
fd(t) = Fs[edge_idx](t) .* Ds[edge_idx](t)
deriv_fd(t) = ForwardDiff.derivative(fd, t)

q1(t) = (2 * λ[1] * E(t)[1] * Ps[edge_idx](t)[1] - pderiv(t)[1] - (λ[1] + μ[1])* Ps[edge_idx](t)[1]) / 
            (Ps[edge_idx](t)[1] - Ps[edge_idx](t)[2])
q2(t) = (2 * λ[2] * E(t)[2] * Ps[edge_idx](t)[2] - pderiv(t)[2] - (λ[2] + μ[2])* Ps[edge_idx](t)[2]) / 
            (Ps[edge_idx](t)[2] - Ps[edge_idx](t)[1])

# fd
q1(t) = (2 * λ[2] * E(t)[2] * fd(t)[2] - deriv_fd(t)[2] - (λ[2] + μ[2])* fd(t)[2]) / 
            (fd(t)[2] - fd(t)[1])

x = collect(range(Fs[edge_idx].t[1], Fs[edge_idx].t[end], length = 50 ))
plot(x, q1.(x))
plot!(x, q2.(x))

using QuadGK
a = Fs[edge_idx].t[1]
b = Fs[edge_idx].t[end]
N1 = quadgk(t -> q1(t), a, b)[1]
N2 = quadgk(t -> q2(t), a, b)[2]




N1 + N2
Nscm[104]


N1 + N2


## PLOT F for branch 104 in cetacea
edge_idx = 104

a = Fs[edge_idx].t[1]
b = Fs[edge_idx].t[end]

x = collect(range(a, b, length =50))

p9 = plot(x, hcat(Fs[edge_idx].(x)...)', ylab = "F(t)", xlab = "time (age)", xflip = true, title = "branch 104")

edge_idx = 1

a = Fs[edge_idx].t[1]
b = Fs[edge_idx].t[end]

x = collect(range(a, b, length =50))

p10 = plot(x, hcat(Fs[edge_idx].(x)...)', ylab = "F(t)", xlab = "time (age)", xflip = true, title = "branch 1")


plot(p9, p10, ylim = (0.0, 0.9))



## Transition probability in a small time span Δt


edge_idx = 104

a = Fs[edge_idx].t[1]
b = Fs[edge_idx].t[end]

ntimes = 30
times = collect(range(a, b, length = ntimes))
Δt = times[2] - times[1]

trans_probs = zeros(ntimes-1, K, K)

for i in 1:(ntimes-1)
    #F0 = [1.0, 0.0]
    #Fₜ₊₁ = 
    #C = ((Ds[edge_idx](times[i]) * ones(K)') .* (I(K) .+ Amatrix(times[i]) * I(K) * Δt))'
    tspan = (times[i], times[i] + Δt)
    p = [λ, μ, η, NaN, K, E]
    F0s = [[1.0, 0.0], [0.0, 1.0]]

    res = zeros(K,K)
    for (j, F0) in enumerate(F0s)
        prob = ODEProblem(Diversification.forward_prob, F0, tspan, p)
        sol = solve(prob)[end]
        res[j,:] = sol
    end

    #Δt .* Amatrix(t) * I(K)
    #C = I(K) .+ Δt .* Amatrix(times[i]) * I(K)
    C = I(K) .+ res

    trans_probs[i,:,:] = C ./ ((C * ones(K)) * ones(K)')
end

trans_probs[:,1,2]


trans_probs[1,:,:] ./ sum(trans_probs[1,:,:], dims = 2)

trans_probs[1,:,:] ./ ((trans_probs[1,:,:] * ones(K))*ones(K)')

trans_probs[1,:,:]

"asd" * 1.0

import Base: *

function *(x::Float64, s::String)
    return(s)
end

*(1.0, "hello")





sum([
    1.0 2.0
    3.0 4.0
], dims = 2)

P_unnormalized = Ds[edge_idx](times[1]) .* ([1.0, 0.0] .+ Δt .* Amatrix(times[1]) * [1.0, 0.0])


trans_probs[1,:,:]

trans_probs[end,:,:] ./ Δt

plot(times[2:end], trans_probs[:,1,2], ylim = (0.0, 0.1))
plot!(times[2:end], trans_probs[:,2,1])

Amatrix(times[1]) * [1.0, 0.0] .* Δt

## Attempt 4

# dN/dt = dP/dt + D2(t)

function my_ode100!(dN, N, p, t)
    λ, μ, E, D, S, K = p

    Q = -I(K) .* η .+ (1 .- I(K)) .* (η/(K-1))
    A = (λ .+ μ .- 2 .* λ .* E(t)) .* I(K) .- Q

    P1 = [1.0, 0.0] .+ A * [1.0, 0.0]
    P1 = P1 ./ sum(P1)

    P2 = [0.0, 1.0] .+ A * [0.0, 1.0]
    P2 = P2 ./ sum(P2)
    
    dN[1] = P1[2] * S(t)[2]
    dN[2] = P2[1] * S(t)[1]
end




edge_idx = 104
a = Fs[edge_idx].t[1]
b = Fs[edge_idx].t[end]
tspan = (a, b)
p = [λ, μ, E, Ds[edge_idx], Ps[edge_idx], K]

N0 = [0.0, 0.0]
prob = ODEProblem(my_ode100!, N0, tspan, p)
solve(prob)
plot(solve(prob))


plot(plot(Fs[104], title = "F"), plot(Ds[104], xflip = true, title = "D"), ylim = (0.0, 0.8))

foo(t) = ForwardDiff.derivative(Fs[edge_idx], t) .- ForwardDiff.derivative(Ds[edge_idx],t)

foo(17)

## compute G(t) as in L+P (2020) sys bio

function G_ode!(dG, G, p, t)
    λ, μ, E, D, S, K = p

    Q = -I(K) .* η .+ (1 .- I(K)) .* (η/(K-1))
    A = diagm(-λ .- μ .+ 2 .* λ .* E(t)) .+ Q

    dG[:,:] = A * G
end

edge_idx = 1
a = Fs[edge_idx].t[end]
b = Fs[edge_idx].t[1]
#tspan = (0.0, maximum(data.branching_times))
tspan = (a, b)
p = [λ, μ, E, Ds[edge_idx], Ps[edge_idx], K]

G0 = zeros(K,K) + I(K)
prob = ODEProblem(G_ode!, G0, tspan, p)
solG = solve(prob, Tsit5(), dt = 0.01)
plot(solG)

## F(t+s) = G(t+s) * F(t)
X0 = solG(a) \ Ds[1](a)

solG(a+2.5) * X0
Ds[1](a+2.5)

(solG(a+2.5) * X0) ./ Ds[1](a+2.5)

solG.t

Fs[1](30)



edge_idx = 104
a = Ds[edge_idx].t[1]
b = Ds[edge_idx].t[end]

foobar1(t) = η * Ds[edge_idx](t)[2] / sum(Ds[edge_idx](t))
foobar2(t) = η * Ds[edge_idx](t)[1] / sum(Ds[edge_idx](t))

n1 = quadgk(foobar1, a, b)[1]
n2 = quadgk(foobar2, a, b)[1]

Nscm[104]

function my_another_ode(dN, N, p, t)
    D, F, K, η = p

    denom = sum(F(t) .* D(t))    

    dN[1] = (η /(K-1)) * (F(t)[2] * D(t)[2]) / denom
    dN[2] = (η /(K-1)) * (F(t)[1] * D(t)[1]) / denom
end


edge_idx = 1
a = Ds[edge_idx].t[1]
b = Ds[edge_idx].t[end]
tspan = (a, b)
p = [Ds[edge_idx], Fs[edge_idx], K, η]
N0 = [0.0, 0.0]

prob1 = ODEProblem(my_another_ode, N0, tspan, p)
solve(prob1)[end]

function normalize_rows(A)
    norm = sum(A, dims = 2)
    res = A ./ (norm * ones(size(A)[1])')
    return(res)
end

normalize_rows(I(K) .+ Δt .* Amatrix(10) * I(K))

function onehot(i, K)
    x = zeros(K)
    x[i] = 1.0
    return(x)
end

Pmatrix(t, Δt) = [
    (onehot(1, K) .+ Δt .* Amatrix(t) * onehot(1, K))[1] (onehot(1, K) .+ Δt .* Amatrix(t) * onehot(1, K))[2]
    (onehot(2, K) .+ Δt .* Amatrix(t) * onehot(2, K))[1] (onehot(2, K) .+ Δt .* Amatrix(t) * onehot(2, K))[2]
]

Pmatrix(10.0, Δt)
## this is equivalent to above
Pm = I(K) .+ Δt .* Amatrix(10.0) * I(K)

0.00072439 / (-Δt)



factorize(Pm)

svd(Pm)


names(df)

ntimes = 1000

edge_idx = 104
a = Ds[edge_idx].t[end]
b = Ds[edge_idx].t[1]
times = collect(range(a, b, length = ntimes))
Δt = times[2] - times[1]

ys = zeros(K, ntimes-1)
for i in 1:(ntimes-1)
    n = abs.(Ps[edge_idx](times[i]) .- Ps[edge_idx](times[i+1]))
    ys[:,i] = n
end

Estatechange = zeros(ntimes-1)
for i in 1:(ntimes-1)
    #n = abs.(Ps[edge_idx](times[i]) .- Ps[edge_idx](times[i+1]))
    Pm = normalize_rows(I(K) .+ Δt .* Amatrix(times[i]) * I(K))
    #fnormalized(t) = Fs[edge_idx]
    current_val = Fs[edge_idx](times[i]) ./ sum(Fs[edge_idx](times[i]))
    Estatechange[i] = ((current_val * ones(K)') .* Pm[1,2])[1,2]
end
Estatechange |> sum

sum(ys, dims = 2)


Nscm[1]






tspan = (a, b)
p = [Ds[edge_idx], Fs[edge_idx], K, η]
N0 = [0.0, 0.0]

function trans_prob8(model, t, Δt)
    K = length(model.λ)

    A = Amatrix(t)

    P_unnorm = I(2) .- Δt .* Amatrix(t) 
    rsum = sum(I(K) .- Δt .* Amatrix(t), dims = 2)
    P = P_unnorm ./ rsum
    return(P)
end

nshifts = zeros(nbranches)
@showprogress for edge_idx in 1:nbranches
    a = Fs[edge_idx].t[1]
    b = Fs[edge_idx].t[end]

    ntimes = 1000
    times = collect(range(a, b, length = ntimes))
    Δt = times[2] - times[1]

    nshift = 0.0
    for i in 2:ntimes
        P = trans_prob8(model, times[i], Δt)
        P[1,1] = 0.0
        P[2,2] = 0.0

        state_prob = Ss[edge_idx](times[i])
        #state_prob = Fs[edge_idx](times[i]) ./ sum(Fs[edge_idx](times[i]))
        #nshift += state_prob[1] .* P[1,2]
        #nshift += state_prob[2] .* P[2,1]
        nshift += sum(P * state_prob)
    end
    nshifts[edge_idx] = nshift
end

start1 = Ds[30](Ds[30].t[end]) .* Ds[1](Ds[1].t[end]) .* λ

Fs[30](Fs[30].t[1]) ./ sum(Fs[30](Fs[30].t[1]))
Ss[30](Fs[30].t[1])

(Ds[1](Ds[1].t[end]) .* Ds[30](Ds[30].t[end]) .* λ) ./ sum(Ds[1](Ds[1].t[end]) .* Ds[30](Ds[30].t[end]) .* λ)

start1 ./ sum(start1)
Ss[30](Ds[30].t[end])

nshifts
plot(nshifts, Nscm, linetype = :scatter)
plot(nshifts, data.branch_lengths, linetype = :scatter)
plot!([nshifts[104]],[ data.branch_lengths[104]], color = "red", lab = "", linetype = :scatter)
annotate!(nshifts, data.branch_lengths, 1:nbranches)

plot(data.branch_lengths, nshifts ./ data.branch_lengths, linetype = :scatter)

histogram(nshifts ./ data.branch_lengths)



phy = Dict("edge" => data.edges,
      "tip.label" => data.tiplab,
      "Nnode" => length(data.tiplab)-1,
     "edge.length" => data.branch_lengths)

## cd /home/bkopper/projects/BDS_deterministic_map

node = data.edges[:,2]
@rput nshifts
@rput phy
@rput node

R"""
library(ape)
library(ggtree)
library(tidytree)
library(ggplot2)
library(dplyr)
"""

R"""
class(phy) <- "phylo"
th <- max(node.depth.edgelength(phy))

df1 <- tibble("node" = node,
            "nshifts" = nshifts)
df1 <- rbind(df1, 
            tibble("node" = length(phy$tip.label),
                   "nshifts" = 0))
x <- as_tibble(phy)

phydf <- merge(x, df1, by = "node")
td_phy <- as.treedata(phydf)

p1a <- ggtree(td_phy, aes(color = `nshifts`)) +
    geom_tiplab(size = 8) +
    theme(legend.position = c(0.2, 0.8)) +
    xlim(c(0.0, th + 10)) 
ggsave("figures/cetacea_2state.pdf", p1a)
""";



