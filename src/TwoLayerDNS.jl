"""
    struct TwoLayerDNS
Container for all the elements of a `TwoLayerDNS`.
"""
struct TwoLayerDNS{NHM <: NonhydrostaticModel, CPF <: ContinuousProfileFunction,
                   TLIC <: TwoLayerInitialConditions, SP <: SalinityPerturbation}
    "An [Oceananigans.jl `NonhydrostaticModel`](https://clima.github.io/OceananigansDocumentation/dev/appendix/library/#Oceananigans.Models.NonhydrostaticModels.NonhydrostaticModel-Tuple{})"
    model :: NHM
    "Continuous profile function"
    profile_function :: CPF
    "The two layer initial conditions"
    initial_conditions :: TLIC
    "Perturbation to the salinity"
    salinity_perturbation :: SP
end
function Base.show(io::IO, tldns::TwoLayerDNS)
    println(io, "TwoLayerDirectNumericalSimulation")
    println(io, " ┣━━━━━━━━━━━━━━━━━ model: $(summary(tldns.model))")
    println(io, " ┣━━━━━━ profile_function: $(typeof(tldns.profile_function))")
    println(io, " ┣━━━━ initial_conditions: $(typeof(tldns.initial_conditions))")
    print(io,   " ┗━ salinity_perturbation: $(typeof(tldns.salinity_perturbation))")
end
