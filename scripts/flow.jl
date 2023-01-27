# flow

function A(λs, μs, η, E, K)
    m = Matrix{Float64}(undef, K, K)
    denom = η/(K-1)

    for i in 1:K, j in 1:K
        if j != i
            m[i,j] = denom
        end
    end

    #for i in 1:K
    #    m[i,i] = - (λ[i] + μ[i] + η) + 2 * λ[i] * E[i]
    #end
    m = m .+ Diagonal(- (λs .+ μs .+ η) .+ 2 .* λ .* Et)

    return m
end

Et = abs.(rand(100))
λs = abs.(rand(100))
μs = abs.(rand(100))
K = 100
η = 0.1

@benchmark A(λs, μs, 0.1, Et, 100)

@benchmark Diagonal(- (λs .+ μs .+ η) .+ 2 .* λs .* Et)
d = Diagonal(- (λs .+ μs .+ η) .+ 2 .* λs .* Et)
#Matrix{Float64}(η / (K-1), K, K)
#@benchmark ones(K,K) .* (η / (K-1))

At = A(λs, μs, 0.1, Et, 100)

v = abs.(rand(100))

@btime At \ v

function odefoo!(dG, G, p, t)
    λ, μ, η, E, K = p
    At = A(λ, μ, η, E, K)

    dG[:] = At * G
end

tspan = (0.0, 1.5)
u0 = Matrix{Float64}(I, K, K)
p = λs, μs, η, Et, K

prob = ODEProblem(odefoo!, u0, tspan, p)

Gt = solve(prob)

Gt[1]

Gt[end]

plot(Gt)


η


Q = [
    -3*η η η η
    η -3*η η η
    η η -3*η η
    η η η -3*η
]

exp(Q .* 0.5)





13.206 * 462 / 1000


