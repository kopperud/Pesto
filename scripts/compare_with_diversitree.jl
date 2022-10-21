using Diversification
import ProgressMeter

#############################
##
##   Load the data files
##
###############################
treefile = "data/bears.tre"
datafile = "data/bears.csv"

phy = readtree(treefile)
ρ = 1.0
data = make_SSEdata(phy, datafile, ρ; include_traits = true)

##############################
##
##   Set up the model
##
###############################

λ = [0.1, 0.5]
μ = [0.05, 0.3]
η = 0.05

model = SSEconstant(λ, μ, η)

## Diversitree
println("Calculating ancestral state marginal probabilities for a state-dependent speciation extinction model (SSE). The SSE model has two states, with unequal diversification rates, and we use equal rates for the state transitions (η = q₀₁ = q₁₀). The rates themselves are time-homogeneous. We are using the bears phylogeny (ntaxa = 8) with fake binary tip data.")
println()
println("Diversitree-BISSE implementation:\n")
asr = anc_state_prob_bisse(treefile, datafile, model)
display(asr)
println()

## New approach
println("New approach, forward-backward pass:\n")
#Ds, Fs = backwards_forwards_pass(model, data; verbose = true);
logL = logL_root(model, data)
D_ends, Ds, sf, E = postorder(model, data; verbose = false);
Fs, F_ends = preorder(model, data, E, D_ends; verbose = false);
res = calculate_tree_rates(data, model, Ds, Fs; verbose = true);

ASP = ancestral_state_probs(data, model, D_ends, F_ends);

println("logL = ", logL)
display(ASP)

p = plot(xlab = "Ancestral state probability\n(present approach)", ylab = "Ancestral state probability (diversitree)", legend = :topleft)
plot!(p, [0.0, 1.0], [0.0, 1.0], lab = "One-to-one", color = "black")
for i in 1:2
    scatter!(p, ASP[:,i], asr[:,i], lab = "State "*string(i))
end

#hcat([res["Ps"][i](Fs[i].t[1]) for i in 1:14]...)'


savefig(p, "figures/compare_with_diversitree.pdf")
