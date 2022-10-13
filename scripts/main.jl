using CSV
using RCall
using DifferentialEquations

include("datatypes.jl")
include("utils.jl")
include("diversitree_bisse.jl")
include("postorder.jl")
include("preorder.jl")
include("ODE.jl")
include("logLroot.jl")

λ = [2.0, 1.0]
μ = [0.5, 0.1]
η = 0.1
ρ = 1.0
k = length(λ)

model = SSEconstant(λ, μ, η)

treefile = "data/bears.tre"
datafile = "data/bears.csv"

println("Diversitree-BISSE model:\n")
asr = anc_state_prob_bisse(treefile, datafile, model)
println("Ancestral state marginal probabilities:")
display(asr)

println("\t\t")
println("Backwards-forwards method:")
phy = readtree(treefile)
data = make_SSEdata(phy, datafile, ρ)
D_ends, sf, E = postorder(model, data)
logL = logL_root(model, data)
println("logL: \t", logL)

ASP = preorder(model, data, E, D_ends)
display(ASP)
println()
