## η versus E[N] across the tree

## How much does the number of shifts depend on the η parameter?

using Distributions
using Diversification

phy = readtree(Diversification.path("primates.tre"))
ρ = 0.67
data = make_SSEdata2(phy, ρ)

λml, μml = estimate_constant_bdp(data)

H = 0.587
n = 6
dλ = LogNormal(log(λml), H)
dμ = LogNormal(log(µml), H)

λquantiles = Diversification.make_quantiles2(dλ, n)
λquantiles = Diversification.make_quantiles(dλ, n)
µquantiles = make_quantiles(dμ, n)

λ, μ = allpairwise(λquantiles, µquantiles)

#ηs = [0.0001, 0.001, 0.01, 0.1, 1.0, 10.0]
ηs = Diversification.lrange(0.001, 10.0; length = 20)

models = [SSEconstant(λ, μ, η) for η in ηs]
logLs = [logL_root(model, data) for model in models]
plot(ηs, logLs, linetype = [:scatter, :line], xscale = :log)



nshifts = zeros(length(data.branch_lengths), length(models))
@showprogress for (i, model) in enumerate(models)
    Ds, Fs = Diversification.backwards_forwards_pass(model, data);
    Ss = Diversification.ancestral_state_probabilities(data, model, Ds, Fs)
    E = extinction_probability(model, data)
    nshift = compute_nshifts(model, data, Ds, Ss; ntimes = 500, ape_order = false)
    nshifts[:,i] = nshift
end

ηtl = ηs .* sum(data.branch_lengths)
sumshifts = sum(nshifts, dims = 1)'
p2 = plot(ηtl, sumshifts, linetype = [:scatter], 
        xlab = "η * treelength", ylab = "Sum of all shifts", label = "E[N|η]",
        title = "Primates tree")
maxdot = maximum([maximum(ηtl), maximum(sumshifts)])
plot!(p2, [1, maxdot], [1, maxdot], label = "one-to-one")
plot!(p2, scale = :log, legend = :topleft)

sumshifts = sum(nshifts["lognormal"][:,:], dims = 1)
p3 = plot((ns .+1) .^2 , sumshifts', linetype = [:scatter, :line], 
        xlab = "State space size (K)", ylab = "Sum of all shifts\n(η = 0.01)",
        title = "Lognormal discrete")
plot!(p3, [(ns[1]+1)^2, (ns[end]+1)^2], [ηtl, ηtl], label = "η * treelength")

plot(p2, p3, ylim = (8, 20))


x, y = parameterizations["lognormal"](8)
scatter(x, y)





