using Distributions

struct SSEconstant <: ContinuousUnivariateDistribution
    λ
    μ
    η
end

struct SSEdata
    state_space
    trait_data
    edges
    tiplab
    node_depth
    ρ
    branch_lengths
    branching_times
    po
end

struct BDconstant <: ContinuousUnivariateDistribution
    λ
    μ
end
