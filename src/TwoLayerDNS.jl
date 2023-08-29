# abstract type TwoLayerDNS end

struct TwoLayerDNS
    "An [Oceananigans.jl `NonhydrostaticModel`](https://clima.github.io/OceananigansDocumentation/dev/appendix/library/#Oceananigans.Models.NonhydrostaticModels.NonhydrostaticModel-Tuple{})"
    model :: NonhydrostaticModel
    "Continuous profile function"
    profile_function <: ContinuousProfileFunction
    "The two layer initial conditions"
    initial_conditions <: TwoLayerInitialConditions
end
