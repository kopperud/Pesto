using Distributions
using StatsPlots
using Diversification
using Optim
using BenchmarkTools
using ProgressMeter

##############################
##
##   Load the data files
##
###############################
treefile = "data/primates.tre"
phy = readtree(treefile)
num_total_species = 367
ρ = length(phy[:tip_label]) / num_total_species
data = make_SSEdata2(phy, ρ)

##############################
##
##   Set up the model
##
###############################

H = 0.587405

d1 = LogNormal(log(0.29), 3 * H)
d2 = LogNormal(log(0.20), 3 * H)
η = 0.1

n = 10
speciation = make_quantiles(d1, n)
extinction = make_quantiles(d2, n)

k = n ^2
λ = zeros(k)
μ = zeros(k)

for (i, (sp, ex)) in enumerate(Iterators.product(speciation, extinction))
    λ[i] = sp
    μ[i] = ex
end
model = SSEconstant(λ, μ, η)
pscatter = scatter(model.λ, model.μ, scale = :log10, xlab = "Speciation rate (λ)", ylab = "Extinction rate (μ)", label ="")
savefig(pscatter, "figures/rates_scatter.svg")

##############################
##
## Calculate the backwards-forwards pass equations
##
###############################

res = birth_death_shift(model, data)


E = extinction_probability(model, data);

Ds, Fs = backwards_forwards_pass(model, data; verbose = true);
Ps = ancestral_state_probabilities(data, model, Ds, Fs);

Fs[1].t


Ps[1](51.88)


@benchmark postorder(model, data, E)
@benchmark postorder_nosave(model, data, E)

logL_root(model, data)
import ForwardDiff

function foo(x, data; n = 10)
    λmean, μmean, η = x

    H = 0.587405

    d1 = LogNormal(log(λmean), 3 * H)
    d2 = LogNormal(log(μmean), 3 * H)

    #n = 10
    speciation = Diversification.make_quantiles2(d1, n)
    extinction = Diversification.make_quantiles2(d2, n)

    k = n^2
    λ = zeros(typeof(λmean), k)
    μ = zeros(typeof(μmean), k)

    for (i, (sp, ex)) in enumerate(Iterators.product(speciation, extinction))
        λ[i] = sp
        μ[i] = ex
    end

    model = SSEconstant(λ, μ, η)
    logL = logL_root(model, data)
    return(logL)
end

x = [0.11, 0.09, 0.05]
foo(x, data)

## fixed η
λs = range(0.009, 0.1, length = 100)
μs = range(0.007, 1.5, length = 100)

z = zeros(length(λs), length(μs))
prog = ProgressMeter.Progress(length(λs)*length(μs), "Calculating surface z");
for (i, λ) in enumerate(λs)
    for (j, μ) in enumerate(μs)
        z[i,j] = foo([λ, μ, 0.1], data)
        ProgressMeter.next!(prog)
    end
end

p1 = surface(λs, μs, z', xlab = "λ", ylab = "μ", zlab = "logL")
savefig(p1, "figures/sse_logLsurface_etafixed.pdf")

## fixed μ
λs = range(0.01, 0.2, length = 100)
ηs = range(0.007, 0.6, length = 100)

z = zeros(length(λs), length(ηs))
prog = ProgressMeter.Progress(length(λs)*length(ηs), "Calculating surface z");
for (i, λ) in enumerate(λs)
    for (j, η) in enumerate(ηs)
        z[i,j] = foo([λ, 0.3, η], data)
        ProgressMeter.next!(prog)
    end
end
z[z .< -750] .= -750
p2 = surface(λs, ηs, z', xlab = "λ", ylab = "η", zlab = "logL")
savefig(p2, "figures/sse_logLsurface_mufixed.pdf")

## fixed λ
μs = range(0.009, 1.0, length = 100)
ηs = range(0.001, 0.1, length = 100)

z = zeros(length(μs), length(ηs))
prog = ProgressMeter.Progress(length(μs)*length(ηs), "Calculating surface z");
for (i, μ) in enumerate(μs)
    for (j, η) in enumerate(ηs)
        z[i,j] = foo([0.1, μ, η], data)
        ProgressMeter.next!(prog)
    end
end
p3 = surface(μs, ηs, z', xlab = "μ", ylab = "η", zlab = "logL", camera = (30, 10))
savefig(p3, "figures/sse_logLsurface_lambdafixed.pdf")

logL_root(model, data)

p2 = surface(λs, μs, z', xlab = "λ", ylab = "μ", zlab = "logL")
scatter3d!(p2, repeat(λs, length(μs)), sort(repeat(μs, length(λs))), vec(z), 
            alpha = 0.3, label = "")

savefig(p1, "figures/sse_logLsurface.pdf")


### Plot grid search as a function of increasing category size

## fixed η
λs = range(0.009, 0.1, length = 24)
μs = range(0.007, 1.5, length = 24)
ns = [32, 16, 12, 8, 6, 4, 3, 2]

prog = ProgressMeter.Progress(length(λs)*length(ns), "Calculating surface z");
z1 = zeros(length(λs), length(μs), length(ns))
jj = Threads.Atomic{Int}(0)
l = Threads.SpinLock()

Threads.@threads for (i, λ) in collect(enumerate(λs))
    for (k, n) in enumerate(ns)
        for (j, μ) in enumerate(μs)
            z1[i,j,k] = foo([λ, μ, 0.1], data; n = n)
        end

        Threads.atomic_add!(jj, 1)
        Threads.lock(l)
        ProgressMeter.update!(prog, jj[])
        Threads.unlock(l) 
    end
end

ps = []
for (k, n) in enumerate(ns)
    z = z1[:,:,k]
    zplot = surface(λs, μs, z', xlab = "λ", ylab = "μ", 
                    zlab = "logL", title = "n=$n", 
                    camera = (30, 10))
    append!(ps, [zplot])
end
pp = plot(ps..., colorbar = :none)
savefig(pp, "figures/more_categories.svg")

ps2 = []
for (k, n) in enumerate(ns)
    z = z1[:,:,k]'
    zc = deepcopy(z)
    mx = maximum(zc)
    zc[zc .< mx - 30] .= mx - 30
    p = StatsPlots.contour(λs, μs, zc, zlim = (maximum(zc)-30, maximum(zc)), 
                           title = "n=$n", levels = 12)
    append!(ps2, [p])
end
pp2 = plot(ps2..., colorbar = :none)
savefig(pp2, "figures/more_categories_contour.svg")


@elapsed foo([0.1, 0.05, 0.1], data; n = 32)

## 
ancestors = Diversification.make_ancestors(data)

ancestor_node = Dict(val => key for (key, val) in eachrow(data.edges))

## Try a random sample from the multivariate lognormal distribution
σ = 0.15
Σ = [0.4 σ
     σ 0.3]
mv = Distributions.MultivariateNormal([2.2, 1.6], Σ)

mvx = rand(mv, 100); scatter(exp.(mvx[1,:]), exp.(mvx[2,:]))

surface(0.0001:0.01:1.8,
        0.0001:0.01:1.8,
        (x1, x2) -> exp.(logpdf(mv, [exp(x1), exp(x2)])), 
        scale = :identity)


logLs = zeros(100)
@showprogress for i in 1:100
    η = 0.1
    mvx = rand(mv, 100)
    model = SSEconstant(exp.(mvx[1,:]), exp.(mvx[2,:]), η)
    logLs[i] = logL_root(model, data)
end

StatsPlots.histogram(logLs, bins = 50)



mv = Distributions.MultivariateNormal([0.1, 0.05], Σ)
results = []
@showprogress for i in 1:100
    mvx = rand(mv, 20)
    model = SSEconstant(exp.(mvx[1,:]), exp.(mvx[2,:]), η)

    res = birth_death_shift(model, data)
    append!(results, [res])
end



p = plot()
for i in 1:100
    violin!(p, [i], [x["lambda"][i] for x in results], bins = 10, label = "", color = "black")
end
p




## Lower bounds of the horseshoe probability density
K = 1/((π^3)^0.5)
phorseshoe1(x) = (K/2) * log(1 + 4 / x^2)
phorseshoe2(x) = K * log(1 + 2 / x^2)


x = -3:0.005:3
plot(x, phorseshoe1.(x), label = "Horseshoe (lower bound)", ylab = "f(x)", xlab = "x")
plot!(x, phorseshoe2.(x), label = "Horseshoe (upper bound)")
plot!(x, pdf.(Normal(0.0, 1.0), x), label = "Normal")
plot!(x, pdf.(Cauchy(0.0, 1.0), x), label = "Cauchy")


Distributions.MultivariateNormal


plot(p1, p2, layout =  (2,1))


scatter3d!(p, λs, μs, z)
p = surface(λs, μs, (λ, μ) -> foo([λ, μ, 0.1], data), xlab = "λ", ylab = "μ")
scatter!(p, λs, μs, (λ, μ) -> foo([λ, μ, 0.1], data), xlab = "λ", ylab = "μ"))


function estimate_sse_hyperparam(data::SSEdata; xinit = [0.11, 0.09, 0.05], lower = [0.0001, 0.0001, 0.0001], upper = [20.0, 20.0, 20.0])
    ρ = data.ρ

    ## ML estimates of parameters
    f(x) = -foo(x, data) ## function to minimize

    inner_optimizer = Optim.GradientDescent()
    optres = Optim.optimize(f, lower, upper, xinit, Fminbox(inner_optimizer))

    λml, μml, ηml = optres.minimizer
    #return(λml, μml, ηml)
    return(optres)
end

@elapsed estimate_sse_hyperparam(data)

using ForwardDiff

ForwardDiff.gradient(foo, [0.2, 0.1, 0.1])

