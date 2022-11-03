using Diversification, DifferentialEquations, Distributions


H = 0.587405
d1 = LogNormal(log(0.08251888), 3 * H)
d2 = LogNormal(log(0.08251888/2.0), 3 * H)

n = 6
speciation = make_quantiles(d1, n)
extinction = make_quantiles(d2, n)
η = 0.1

k = n ^2
λ = zeros(k)
μ = zeros(k)

for (i, (sp, ex)) in enumerate(Iterators.product(speciation, extinction))
    λ[i] = sp
    μ[i] = ex
end

model = SSEconstant(λ, μ, η)

i_not_js = [setdiff(1:k, i) for i in 1:k]
pE = [model.λ, model.μ, model.η, i_not_js, k]
tspan = (0.0, 101.35870000000006)

E0 = repeat([0.0], k)
pr = DifferentialEquations.ODEProblem(Diversification.extinction_prob, E0, tspan, pE);
E = DifferentialEquations.solve(pr, alg)

pr_old = DifferentialEquations.ODEProblem(Diversification.extinction_prob_old, E0, tspan, pE);
E_old = DifferentialEquations.solve(pr_old, alg);

pD = [model.λ, model.μ, model.η, i_not_js, k, E]
D0 = repeat([1.0], k)
 
prob = DifferentialEquations.ODEProblem(Diversification.backward_prob, D0, tspan, pD);
D = DifferentialEquations.solve(prob, alg)


