using ParameterizedFunctions, OrdinaryDiffEq,
      ODEInterfaceDiffEq, Plots, Sundials, SciPyDiffEq, deSolveDiffEq
#using MATLABDiffEq
using DiffEqDevTools
using LinearAlgebra, StaticArrays

#f = @ode_def_bare bisse_extinction begin
#    dE[:] .= μ .- (λ .+ μ .+ η) .* E .+ λ .* E.^2 .+ (η/(K-1)) .* (sum(E) .- E)
#end λ μ η
f = @ode_def_bare bisse_extinction begin
    #λ1, λ2, μ1, μ2, η, K = p
#    dE[:] .= μ .- (λ .+ μ .+ η) .* E .+ λ .* E.^2 .+ (η/(K-1)) .* (sum(E) .- E)
    dE1 = μ1 - (λ1 + μ1 + η) * E1 + λ1 * E1^2 + η * E2
    dE2 = μ2 - (λ2 + μ2 + η) * E2 + λ2 * E2^2 + η * E1
end λ1 λ2 μ1 μ2 η K
λ = [0.1, 0.2]
μ = [0.05, 0.15]
η = 0.05
K = 2
p = [λ..., μ..., η, K]

tspan = (0.0,10.0)
u0 = [0.0, 0.0]
prob = ODEProblem(f,u0,tspan,p)
#staticprob = ODEProblem{false}(f,u0,tspan,p)
staticprob = ODEProblem{false}(f,SVector{2}(u0),tspan,SVector{6}(p))

sol = solve(prob,Vern7(),abstol=1/10^14,reltol=1/10^14)
test_sol = TestSolution(sol)

## solve using Euler as well
eulers = []
euler_times = []

dts = [0.1^n for n in 1:6]
for dt in dts
    prob = ODEProblem(f,u0,tspan,p, dt = dt)
    ts = []
    for i in 1:3
        euler_time = @elapsed solve(prob, Euler())
        append!(ts, euler_time)
    end
    euler_time = minimum(ts)
    euler_sol = solve(prob, Euler())
    append!(eulers, [euler_sol])
    append!(euler_times, euler_time)
end

#euler_errors = [sum(abs.(e.u[end] - sol.u[end])) for e in eulers]
euler_errors = [LinearAlgebra.norm(Statistics.mean(e.u .- sol.(e.t))) for e in eulers]
#euler_errors = [LinearAlgebra.norm(Statistics.mean(e.(sol.t) .- sol.u)) for e in eulers]




setups = [
          Dict(:alg=>DP5())
          Dict(:alg=>Tsit5())
          Dict(:alg=>Vern7())
          Dict(:alg=>RK4())
          Dict(:prob_choice => 2, :alg=>DP5())
          Dict(:prob_choice => 2, :alg=>Tsit5())
          Dict(:prob_choice => 2, :alg=>Vern7())
          Dict(:prob_choice => 2, :alg=>RK4())
          Dict(:alg=>dopri5())
          Dict(:alg=>ODEInterfaceDiffEq.dop853())
          Dict(:alg=>ODEInterfaceDiffEq.radau5())
          Dict(:alg=>ODEInterfaceDiffEq.rodas())
#          Dict(:alg=>MATLABDiffEq.ode45())
#          Dict(:alg=>MATLABDiffEq.ode113())
          Dict(:alg=>SciPyDiffEq.RK45())
          Dict(:alg=>SciPyDiffEq.LSODA())
          Dict(:alg=>SciPyDiffEq.odeint())
          Dict(:alg=>deSolveDiffEq.lsoda())
          Dict(:alg=>deSolveDiffEq.ode45())
          Dict(:alg=>deSolveDiffEq.radau())
          Dict(:alg=>CVODE_Adams())
          Dict(:alg=>CVODE_BDF())
  ]

labels = [
  "Julia: DP5"
  "Julia: Tsit5"
  "Julia: Vern7"
  "Julia: RK4"
  "Julia: DP5 Static"
  "Julia: Tsit5 Static"
  "Julia: Vern7 Static"
  "Julia: RK4 Static"
  "Hairer: dopri5"
  "Hairer: dop853"
  "Hairer: radau5"
  "Hairer: rodas"
#  "MATLAB: ode45"
#  "MATLAB: ode113"
  "SciPy: RK45"
  "SciPy: LSODA"
  "SciPy: odeint"
  "deSolve: lsoda"
  "deSolve: ode45"
  "deSolve: radau"
  "Sundials: Adams"
  "Sundials: BDF"
  ]

abstols = 1.0 ./ 10.0 .^ (6:13)
reltols = 1.0 ./ 10.0 .^ (3:10)
wp = WorkPrecisionSet([prob,staticprob],abstols,reltols,setups;
                      names = labels,print_names = true,
                      appxsol=[test_sol,test_sol],dense=false,
                      save_everystep=false,numruns=100,maxiters=10000000,
                      timeseries_errors=false,verbose=false)
benchmark_plot = plot(wp,
                      title="E(t). λ = [0.1, 0.2], μ = [0.05, 0.15], η = 0.05, from t₀=0.0 to tₕ = 10.0",
                      legend=:outertopleft,
                      color=permutedims([repeat([:LightGreen],4)...,
                                         repeat([:DarkGreen],4)...,
                                         repeat([:Red], 4)...,
                                         #repeat([:Orange],2)...,
                                         repeat([:Yellow],3)...,
                                         repeat([:Blue],3)...,
                                         repeat([:Purple],2)...]),
                      size = (1400,650),
                      legendfontsize = 15,
         xticks = 10.0 .^ (-12:1:5),
         yticks = (10.0 .^ (-5.5:0.5:2), ["3.2 µs", "10 µs", "32 µs", "100 µs", "316 µs", "1 ms", "32 ms", "100 ms", "316 ms", "1s", "3.2 s", "10 s", "32 s"]),
     bottom_margin = 5Plots.mm)

plot!(benchmark_plot, euler_errors, euler_times, label = "", color = "black")
scatter!(benchmark_plot, euler_errors, euler_times, label = "Euler's method", color = "black")
for i in 1:length(dts)
    annotate!(benchmark_plot, euler_errors[i], euler_times[i], text("dt = 0.1^"*string(i), 10, :right ,:top))
end

#plot!(benchmark_plot, xflip = true)


savefig(benchmark_plot, "figures/benchmark_odesolvers.pdf")
savefig(benchmark_plot, "figures/benchmark_odesolvers.png")


