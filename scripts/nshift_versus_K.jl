## How much does the number of shifts depend on the η parameter?

## How much does the number of shifts depend on the size of the state space?



using Distributions
using Diversification

phy = readtree(Diversification.path("primates.tre"))
ρ = 0.67
data = make_SSEdata2(phy, ρ)

λml, μml = estimate_constant_bdp(data)



ns = [2, 3, 4, 5, 6, 7, 8, 9, 10]

parameterizations = Dict()

parameterizations["linear"] = function(n)
    λs = collect(range(0.1, 0.6; length = n +1))
    μs = collect(range(0.0, 0.5; length = n +1))

    λ, μ = allpairwise(λs, µs)
    return(λ, μ)
end

parameterizations["lognormal"] = function(n)
    H = 0.587
    dλ = LogNormal(log(λml), H)
    dμ = LogNormal(log(µml), H)

    λquantiles = make_quantiles(dλ, n)
    µquantiles = make_quantiles(dμ, n)

    λ, μ = allpairwise(λquantiles, µquantiles)
    return(λ, μ)
end

parameterizations["lognormal2"] = function(n)
    H = 0.587
    dλ = LogNormal(log(λml), H)
    dμ = LogNormal(log(µml), H)

    λquantiles = Diversification.make_quantiles3(dλ, n)
    µquantiles = Diversification.make_quantiles3(dμ, n)

    λ, μ = allpairwise(λquantiles, µquantiles)
    return(λ, μ)
end

models = Dict()
for (key, val) in parameterizations
    ms = []
    for n in ns
        λ, μ = val(n)
        η = 0.01

        model = SSEconstant(λ, μ, η)
        append!(ms, [model])
    end
    models[key] = ms
end

logLs = Dict()
for (key, vals) in models
    logLs[key] = [logL_root(model, data) for model in vals]
end

logLs


plot(ns.^2, logLs["linear"], linetype = [:scatter, :line], xlab = "State space (K)", ylab = "logL")
plot!(ns.^2, logLs["lognormal"], linetype = [:scatter, :line], xlab = "State space (K)", ylab = "logL")


nshifts = Dict()
for (key, ms) in models
    x = zeros(length(data.branch_lengths), length(ns))
    @showprogress for (j, model) in enumerate(ms)
        Ds, Fs = Diversification.backwards_forwards_pass(model, data);
        Ss = Diversification.ancestral_state_probabilities(data, model, Ds, Fs)
        E = extinction_probability(model, data)
        nshift = compute_nshifts(model, data, Ds, Ss; ntimes = 100, ape_order = false)
        x[:,j] = nshift
    end
    nshifts[key] = x
end

ηtl = η * sum(data.branch_lengths)
sumshifts = sum(nshifts["linear"][:,:], dims = 1)
p2 = plot(ns .^2 , sumshifts', linetype = [:scatter, :line], 
        xlab = "State space size (K)", ylab = "Sum of all shifts\n(η = 0.01)",
        title = "linear spacing")
plot!(p2, [ns[1]^2, ns[end]^2], [ηtl, ηtl], label = "η * treelength")

sumshifts = sum(nshifts["lognormal"][:,:], dims = 1)
p3 = plot((ns .+1) .^2 , sumshifts', linetype = [:scatter, :line], 
        xlab = "State space size (K)", ylab = "Sum of all shifts\n(η = 0.01)",
        title = "Lognormal discrete")
plot!(p3, [(ns[1]+1)^2, (ns[end]+1)^2], [ηtl, ηtl], label = "η * treelength")

sumshifts = sum(nshifts["lognormal2"][:,:], dims = 1)
p4 = plot((ns .+1) .^2 , sumshifts', linetype = [:scatter, :line], 
        xlab = "State space size (K)", ylab = "Sum of all shifts\n(η = 0.01)",
        title = "Lognormal discrete2")
plot!(p4, [(ns[1]+1)^2, (ns[end]+1)^2], [ηtl, ηtl], label = "η * treelength")

plot(p2, p3, p4, ylim = (8, 20), size = (700, 700))


x, y = parameterizations["lognormal"](8)
scatter(x, y)





