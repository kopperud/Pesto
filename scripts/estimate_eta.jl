
## estimate eta

using Diversification
using Distributions
using Plots

phy = readtree(Diversification.path("primates.tre")) 
ρ = 0.67  
data = make_SSEdata2(phy, ρ)

## estimate constant-rate extinction and speciation rate
λml, μml = estimate_constant_bdp(data)


H = 0.587405
d1 = LogNormal(log(λml), H)
d2 = LogNormal(log(μml), H)

n = 6
speciation = make_quantiles(d1, n)
extinction = make_quantiles(d2, n)

λ, μ = allpairwise(speciation, extinction);


## estimate η conditional on λ, μ
ηml = optimize_eta(λ, μ, data)

Diversification.optimize_eta(λ, μ, data)
sselp(η,  λ, μ, )



#xs = 0.001:0.01:1.0
xs = range(0.001, 10, length = 100)
#ys = [Diversification.f(y, λ, μ, data) for y in xs]
ys = [Diversification.sselp(y, λ, μ, data) for y in xs]

Diversification.sselp(10.0, λ, μ, data)
Diversification.sselp(ηml, λ, μ, data)

plot(xs, ys, xlab = "η", ylab = "logL", linetype = [:scatter, :line])

import Optim
Optim.optimize()




