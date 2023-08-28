"""
    module TwoLayerDNS
Module containing the setup for a two layer Direct Numerical Simulation. This two layer
model is mainly being used to experiments to investigate the cabbeling instability. The
constants that are exported need not be used, they are just values that I have chosen to use
for experiments, see [CabbelingExperiments](https://github.com/jbisits/CabbelingExperiments).
"""
module TwoLayerDNS

using DirectNumericalCabbelingShenanigans, JLD2, GibbsSeaWater
using DirectNumericalCabbelingShenanigans: simulation_progress
using SpecialFunctions: erf
using Oceanostics: KineticEnergyDissipationRate
import Base: show

export
    StableUpperLayerInitialConditions,
    CabbelingUpperLayerInitialConditions,
    UnstableUpperLayerInitialConditions,
    IsohalineUpperLayerInitialConditions,
    IsothermalUpperLayerInitialConditions,
    TwoLayerInitialConditions,
    StableTwoLayerInitialConditions,
    CabbelingTwoLayerInitialConditions,
    UnstableTwoLayerInitialConditions,
    IsohalineTwoLayerInitialConditions,
    IsothermalTwoLayerInitialConditions,
    HyperbolicTangent, Erf,
    GaussianProfile,
    GaussianBlob,
    RandomPerturbations,
    set_two_layer_initial_conditions!,
    add_velocity_random_noise!,
    DOMAIN_EXTENT,
    HIGH_RESOLUTION,
    SO_DIFFUSIVITIES,
    REFERENCE_DENSITY,
    INTERFACE_LOCATION,
    non_dimensional_numbers,
    compute_density,
    compute_density!

"""
    abstract type UpperLayerInitialConditions
Abstract super type for initial temperature and salinity in the upper layer.
"""
abstract type UpperLayerInitialConditions end
"`show` for `UpperLayerInitialConditions`"
function Base.show(io::IO, ulic::UpperLayerInitialConditions)
    if ulic isa IsothermalUpperLayerInitialConditions
        println(io, "$(typeof(ulic))")
        println(io, " ┣━━ S = $(ulic.S₀ᵘ)")
        print(io,   " ┗━━ T = $(ulic.T)")
    elseif ulic isa IsohalineUpperLayerInitialConditions
        println(io, "$(typeof(ulic))")
        println(io, " ┣━━ S = $(ulic.S)")
        print(io,   " ┗━━ T = $(ulic.T₀ᵘ)")
    else
        println(io, "$(typeof(ulic))")
        println(io, " ┣━━ S = $(ulic.S₀ᵘ)")
        print(io,   " ┗━━ T = $(ulic.T₀ᵘ)")
    end
end
"""
    struct StableUpperLayerInitialConditions
Container for initial salinity and temperature conditions that are stable relative to `S₀ˡ`
and `T₀ˡ`.
"""
struct StableUpperLayerInitialConditions{T} <: UpperLayerInitialConditions
    "Initial salinity in the upper layer"
    S₀ᵘ :: T
    "Initial temperature in the upper layer"
    T₀ᵘ :: T
end
"""
    struct CabbelingUpperLayerInitialConditions
Container for initial salinity and temperature conditions that are unstable to cabbeling
relative to `S₀ˡ` and `T₀ˡ`.
"""
struct CabbelingUpperLayerInitialConditions{T} <: UpperLayerInitialConditions
    "Initial salinity in the upper layer"
    S₀ᵘ :: T
    "Initial temperature in the upper layer"
    T₀ᵘ :: T
end
"""
    struct UnstableUpperLayerInitialConditions
Container for initial salinity and temperature conditions that are unstable relative to `S₀ˡ`
and `T₀ˡ`.
"""
struct UnstableUpperLayerInitialConditions{T} <: UpperLayerInitialConditions
    "Initial salinity in the upper layer"
    S₀ᵘ :: T
    "Initial temperature in the upper layer"
    T₀ᵘ :: T
end
"""
    struct IsohalineUpperLayerInitialConditions
Container for isohaline initial salinity at (`S₀ˡ`) and initial temperature conditions `T₀ˡ`.
"""
struct IsohalineUpperLayerInitialConditions{T} <: UpperLayerInitialConditions
    "Initial (uniform) salinity over the domain"
    S   :: T
    "Initial temperature in the upper layer"
    T₀ᵘ :: T
end
"""
    struct IsohalineUpperLayerInitialConditions
Container for isohaline initial salinity at (`S₀ˡ`) and initial temperature conditions `T₀ˡ`.
"""
struct IsothermalUpperLayerInitialConditions{T} <: UpperLayerInitialConditions
    "Initial salinity in the upper layer"
    S₀ᵘ :: T
    "Initial (uniform) temperature over the domain"
    T   :: T
end
"""
    abstract type TwoLayerInitialConditions
Abstract supertype for two layer model initial conditions. Keyword arguments for the
`TwoLayerInitialConditions` constructor `S₀ˡ` and `T₀ˡ` can be specified for the salinity
and temperature of the lower layer but have default values `S₀ˡ = 34.7`gkg⁻¹ and
`T₀ˡ = 0.5`°C.
"""
abstract type TwoLayerInitialConditions end
"`show` for `TwoLayerInitialConditions`"
function Base.show(io::IO, tlic::TwoLayerInitialConditions)
    println(io, "$(typeof(tlic))")
    println(io, " ┣━━ upper_layer: S = $(tlic.S₀ᵘ), T = $(tlic.T₀ᵘ)")
    print(io,   " ┗━━ lower_layer: S = $(tlic.S₀ˡ), T = $(tlic.T₀ˡ)")
end
"""
    struct StableTwoLayerInitialConditions
Container for initial salinity and temperature conditions that are stable.
"""
struct StableTwoLayerInitialConditions{T} <: TwoLayerInitialConditions
    "Initial salinity in the upper layer"
    S₀ᵘ :: T
    "Initial salinity in the lower layer"
    S₀ˡ :: T
    "Initial difference in salinity between the layers"
    ΔS₀ :: T
    "Initial temperature in the upper layer"
    T₀ᵘ :: T
    "Initial temperature in the lower layer"
    T₀ˡ :: T
    "Initial temperature difference between the layers"
    ΔT₀ :: T
end
TwoLayerInitialConditions(initial_conditions::StableUpperLayerInitialConditions;
                          S₀ˡ = 34.7, T₀ˡ = 0.5) =
    StableTwoLayerInitialConditions(initial_conditions.S₀ᵘ, S₀ˡ,
                                    initial_conditions.S₀ᵘ - S₀ˡ, initial_conditions.T₀ᵘ,
                                    T₀ˡ, initial_conditions.T₀ᵘ -T₀ˡ)
"""
    struct CabbelingTwoLayerInitialConditions
Container for initial salinity and temperature conditions that are unstable to cabbeling.
"""
struct CabbelingTwoLayerInitialConditions{T} <: TwoLayerInitialConditions
    "Initial salinity in the upper layer"
    S₀ᵘ :: T
    "Initial salinity in the lower layer"
    S₀ˡ :: T
    "Initial difference in salinity between the layers"
    ΔS₀ :: T
    "Initial temperature in the upper layer"
    T₀ᵘ :: T
    "Initial temperature in the lower layer"
    T₀ˡ :: T
    "Initial temperature difference between the layers"
    ΔT₀ :: T
end
TwoLayerInitialConditions(initial_conditions::CabbelingUpperLayerInitialConditions;
                          S₀ˡ = 34.7, T₀ˡ = 0.5) =
    CabbelingTwoLayerInitialConditions(initial_conditions.S₀ᵘ, S₀ˡ,
                                       initial_conditions.S₀ᵘ - S₀ˡ, initial_conditions.T₀ᵘ,
                                       T₀ˡ, initial_conditions.T₀ᵘ -T₀ˡ)
"""
    struct UnstableTwoLayerInitialConditions
Container for initial salinity and temperature conditions that are gravitationally unstable.
"""
struct UnstableTwoLayerInitialConditions{T} <: TwoLayerInitialConditions
    "Initial salinity in the upper layer"
    S₀ᵘ :: T
    "Initial salinity in the lower layer"
    S₀ˡ :: T
    "Initial difference in salinity between the layers"
    ΔS₀ :: T
    "Initial temperature in the upper layer"
    T₀ᵘ :: T
    "Initial temperature in the lower layer"
    T₀ˡ :: T
    "Initial temperature difference between the layers"
    ΔT₀ :: T
end
TwoLayerInitialConditions(initial_conditions::UnstableUpperLayerInitialConditions;
                          S₀ˡ = 34.7, T₀ˡ = 0.5) =
    UnstableTwoLayerInitialConditions(initial_conditions.S₀ᵘ, S₀ˡ,
                                      initial_conditions.S₀ᵘ - S₀ˡ, initial_conditions.T₀ᵘ,
                                      T₀ˡ, initial_conditions.T₀ᵘ -T₀ˡ)
"""
    struct IsohalineTwoLayerInitialConditions
Container for initial salinity and temperature conditions where salinity is uniform over
domain.
"""
struct IsohalineTwoLayerInitialConditions{T} <: TwoLayerInitialConditions
    "Initial salinity in the upper layer"
    S₀ᵘ :: T
    "Initial salinity in the lower layer"
    S₀ˡ :: T
    "Initial difference in salinity between the layers"
    ΔS₀ :: T
    "Initial temperature in the upper layer"
    T₀ᵘ :: T
    "Initial temperature in the lower layer"
    T₀ˡ :: T
    "Initial temperature difference between the layers"
    ΔT₀ :: T
end
TwoLayerInitialConditions(initial_conditions::IsohalineUpperLayerInitialConditions;
                          T₀ˡ = 0.5) =
    IsohalineTwoLayerInitialConditions(initial_conditions.S, initial_conditions.S,
                                       0.0, initial_conditions.T₀ᵘ, T₀ˡ,
                                       initial_conditions.T₀ᵘ -T₀ˡ)
                                       """
    struct IsothermalTwoLayerInitialConditions
Container for initial salinity and temperature conditions where temperature is uniform over
the domain.
"""
struct IsothermalTwoLayerInitialConditions{T} <: TwoLayerInitialConditions
    "Initial salinity in the upper layer"
    S₀ᵘ :: T
    "Initial salinity in the lower layer"
    S₀ˡ :: T
    "Initial difference in salinity between the layers"
    ΔS₀ :: T
    "Initial temperature in the upper layer"
    T₀ᵘ :: T
    "Initial temperature in the lower layer"
    T₀ˡ :: T
    "Initial temperature difference between the layers"
    ΔT₀ :: T
end
TwoLayerInitialConditions(initial_conditions::IsothermalUpperLayerInitialConditions;
                          S₀ˡ = 34.7) =
    IsothermalTwoLayerInitialConditions(initial_conditions.S₀ᵘ, S₀ˡ,
                                        initial_conditions.S₀ᵘ - S₀ˡ, initial_conditions.T,
                                        initial_conditions.T, 0.0)
"""
    abstract type ContinuousProfileFunction end
Abstract super type for the continuous function that sets the continuous profile for
temperature and salinity.
"""
abstract type ContinuousProfileFunction end
"`show` for `ContinuousProfileFunction`"
function Base.show(io::IO, cpf::ContinuousProfileFunction)
    if cpf isa HyperbolicTangent
        println(io, "$(typeof(cpf))")
        println(io, " ┣━━━━━━━━━━ interface_location: z = $(cpf.interface_location)m ")
        print(io,   " ┗━━ interface_transition_scale: $(cpf.interface_transition_width)")
    elseif cpf isa Erf
        println(io, "$(typeof(cpf))")
        println(io, " ┣━━ interface_location: z = $(cpf.interface_location)m")
        print(io,   " ┗━━━━━━━━ time_evolved: t = $(cpf.time)s")
    end
end
"""
    struct HyperbolicTangent
Container for a hyperbolic tangent profile. The `interface_transition_width` sets the width
of the transition between the upper and lower layer.
"""
struct HyperbolicTangent{T} <: ContinuousProfileFunction
    "Location of the interface between the two layers."
    interface_location :: T
    "Scale the transition between the upper and lower layer salinity and temperature."
    interface_transition_width :: T
end
"""
    struct Erf
Container for a profile that is an error function. The `time` is the time which to evaluate
`erf_tracer_solution`.
"""
struct Erf{T} <: ContinuousProfileFunction
    "Location of the interface between the two layers."
    interface_location :: T
    "Time at which to evaluate the error function which is solution to 1D evolution of S or T."
    time :: T
end
"""
    abstract type SalinityPerturbation
Abstract super type for the salinity perturbation added to the upper layer.
"""
abstract type SalinityPerturbation end
"`show` for `SalinityPerturbation`"
function Base.show(io::IO, sp::SalinityPerturbation)
    if sp isa GaussianProfile
        println(io, "$(typeof(sp))")
        println(io, " ┣━━ interface_location: z = $(sp.interface_location)m")
        println(io, " ┣━━ centre_of_Gaussian: z = $(sp.μ)m")
        println(io, " ┣━━━ width_of_Gaussian: $(sp.σ)")
        print(io,   " ┗━━━━━━━━━━━━━━━ scale: $(sp.scale)")
    elseif sp isa GaussianBlob
        println(io, "$(typeof(sp))")
        println(io, " ┣━━━━━━ depth_location: z = $(round(sp.depth; digits = 3))m")
        println(io, " ┣━━ centre_of_Gaussian: $(sp.μ)")
        println(io, " ┣━━━ width_of_Gaussian: $(sp.σ)")
        print(io,   " ┗━━━━━━━━━━━━━━━ scale: $(sp.scale)")
    elseif sp isa RandomPerturbations
        println(io, "$(typeof(sp))")
        println(io, " ┣━━ noise_location: z = $(round(sp.depth; digits = 3))m")
        print(io,   " ┗━━━━━ noise_scale: $(sp.scale)")
    end
end
"""
    struct GaussianProfile
Container for a Gaussian profile salinity perturbation in the upper layer.
"""
struct GaussianProfile{T} <: SalinityPerturbation
    "Location of interface betweenn upper and lower layers."
    interface_location :: T
    "Center of Gaussin profile in the upper layer."
    μ :: T
    "With of Gaussian profile in the upper layer."
    σ :: T
    "Scale the Gausssian profile, defaults to 1 i.e. it is a pdf."
    scale :: T
end
GaussianProfile(interface_location, μ, σ; scale = 1.0) =
    GaussianProfile(interface_location, μ, σ, scale)
"""
    struct GaussianBlob
Container for a horizontal Gaussian blob salinity perturbation at `depth` in the upper layer.
"""
struct GaussianBlob{T} <: SalinityPerturbation
    "Depth at which to set the horizontal Gaussian blob of salinity."
    depth :: T
    "Centre of the blob."
    μ :: Vector{T}
    "Width of the blob."
    σ :: T
    "Scale the Gausssian blob, defaults to 1 i.e. it is a pdf."
    scale :: T
end
GaussianBlob(interface_location, μ, σ; scale = 1.0) =
    GaussianBlob(interface_location, μ, σ, scale)
"""
    struct RandomPerturbations
Container for adding `scale`d random noise to the salinity at `depth`.
"""
struct RandomPerturbations{T} <: SalinityPerturbation
    "Depth at which to set random salinity perturbations."
    depth :: T
    "Scale for the random noise."
    scale :: T
end
"""
    const DOMAIN_EXTENT
Domain extent on which the two layer simulations are run.
"""
const DOMAIN_EXTENT = (Lx = 0.1, Ly = 0.1, Lz = 1)
"""
    const HIGH_RESOLUTION
Resolution (high) at which to run the DNS.
"""
const HIGH_RESOLUTION = (Nx = 50, Ny = 50, Nz = 1400)
"""
    const SO_DIFFUSIVITIES
Diffusivity estimates for the Southern Ocean.
"""
const SO_DIFFUSIVITIES = (ν = 1e-6, κ = (S = 1e-9, T = 1e-7))
"""
    const REFERENCE_DENSITY
Reference density for use in the two layer DNS. Calculated using the salinity `S₀ˡ` and
temperature `T₀ˡ` of the lower layer .
"""
const REFERENCE_DENSITY = gsw_rho(S₀ˡ, T₀ˡ, 0)
"""
    const INTERFACE_LOCATION
Location of the interface (in the vertical) between the upper and lower layers.
"""
const INTERFACE_LOCATION = -0.375
"""
    function set_two_layer_initial_conditions(model::Oceananigans.AbstractModel,
                                              initial_conditions::TwoLayerInitialConditions,
                                              profile_function::ContinuousProfileFunction)
    function set_two_layer_initial_conditions(model::Oceananigans.AbstractModel,
                                              initial_conditions::TwoLayerInitialConditions,
                                              profile_function::ContinuousProfileFunction,
                                              salinity_perturbation::SalinityPerturbaion)

Set initial conditions for a two layer model that are smooth depending on the
`profile_funciton` with or without a `salinity_perturbation` in the upper layer.

## Function arguments:

- `model`: to set the initial salinity and temperature in;
- `initial_conditions`: the values for the initial conditions in an appropriate container;
- `profile_function`: the smooth funtion used to set the profile to avoid discontinuities;
- `salinity_perturbation`: perturbation for the salinity in the upper layer to help kick off
instability.
"""
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::Erf)

    κₛ, κₜ = model.closure.κ.S, model.closure.κ.T
    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = erf_tracer_solution(z, S₀, ΔS, κₛ, profile_function)
    initial_T_profile(x, y, z) = erf_tracer_solution(z, T₀, ΔT, κₜ, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::Erf,
                                           salinity_perturbation::GaussianProfile)

    κₛ, κₜ = model.closure.κ.S, model.closure.κ.T
    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = erf_tracer_solution(z, S₀, ΔS, κₛ, profile_function) +
                                 perturb_salinity(z, salinity_perturbation)

    initial_T_profile(x, y, z) = erf_tracer_solution(z, T₀, ΔT, κₜ, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::Erf,
                                           salinity_perturbation::GaussianProfile,
                                           salinity_noise::RandomPerturbations)

    κₛ, κₜ = model.closure.κ.S, model.closure.κ.T
    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = erf_tracer_solution(z, S₀, ΔS, κₛ, profile_function) +
                                 perturb_salinity(z, salinity_perturbation) +
                                 perturb_salinity(z, salinity_noise)

    initial_T_profile(x, y, z) = erf_tracer_solution(z, T₀, ΔT, κₜ, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                          initial_conditions::TwoLayerInitialConditions,
                                          profile_function::Erf,
                                          salinity_perturbation::GaussianBlob)

    κₛ, κₜ = model.closure.κ.S, model.closure.κ.T
    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = erf_tracer_solution(z, S₀, ΔS, κₛ, profile_function) +
                                 perturb_salinity(x, y, z, salinity_perturbation)

    initial_T_profile(x, y, z) = erf_tracer_solution(z, T₀, ΔT, κₜ, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::Erf,
                                           salinity_perturbation::RandomPerturbations)

    κₛ, κₜ = model.closure.κ.S, model.closure.κ.T
    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = erf_tracer_solution(z, S₀, ΔS, κₛ, profile_function) +
                                 perturb_salinity(z, salinity_perturbation)

    initial_T_profile(x, y, z) = erf_tracer_solution(z, T₀, ΔT, κₜ, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::HyperbolicTangent)

    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = tanh_initial_condition(z, S₀, ΔS, profile_function)
    initial_T_profile(x, y, z) = tanh_initial_condition(z, T₀, ΔT, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::HyperbolicTangent,
                                           salinity_perturbation::GaussianProfile)

    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = tanh_initial_condition(z, S₀, ΔS, profile_function) +
                                 perturb_salinity(z, salinity_perturbation)

    initial_T_profile(x, y, z) = tanh_initial_condition(z, T₀, ΔT, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::HyperbolicTangent,
                                           salinity_perturbation::GaussianProfile,
                                           salinity_noise::RandomPerturbations)

    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = tanh_initial_condition(z, S₀, ΔS, profile_function) +
                                 perturb_salinity(z, salinity_perturbation) +
                                 perturb_salinity(z, salinity_noise)

    initial_T_profile(x, y, z) = tanh_initial_condition(z, T₀, ΔT, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::HyperbolicTangent,
                                           salinity_perturbation::GaussianBlob)

    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = tanh_initial_condition(z, S₀, ΔS, profile_function) +
                                 perturb_salinity(x, y, z, salinity_perturbation)

    initial_T_profile(x, y, z) = tanh_initial_condition(z, T₀, ΔT, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::HyperbolicTangent,
                                           salinity_perturbation::RandomPerturbations)

    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = tanh_initial_condition(z, S₀, ΔS, profile_function) +
                                 perturb_salinity(z, salinity_perturbation)

    initial_T_profile(x, y, z) = tanh_initial_condition(z, T₀, ΔT, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
"""
    function add_velocity_random_noise!(model::Oceananigans.AbstractModel, noise_magnitude::Number,
                                        location::Number; horizontal = false)
Add standard normally distributed random noise, scaled by `noise_magnitude`, to the
horizontal velocity fields at the interface of the upper and lower layers. If
`location` for where noise should be seeded is not provided then the noise is added
everywhere in the domain. To only add horizontal random noise (i.e. in the `u` and `v`
velocity fields) set `true` for the `horizontal` keyword argument.
"""
function add_velocity_random_noise!(model::Oceananigans.AbstractModel,
                                    noise_magnitude::Number, location::Number;
                                    horizontal = false)

    add_noise(x, y, z) = round(z; digits = 4) == location ? noise_magnitude * randn() : 0

    horizontal == true ? set!(model, u = add_noise, v = add_noise) :
                         set!(model, u = add_noise, v = add_noise, w = add_noise)

    return nothing

end
function add_velocity_random_noise!(model::Oceananigans.AbstractModel,
                                    noise_magnitude::Number; horizontal = false)

    add_noise(x, y, z) = noise_magnitude * randn()

    horizontal == true ? set!(model, u = add_noise, v = add_noise) :
                         set!(model, u = add_noise, v = add_noise, w = add_noise)

    return nothing

end

"""
    function erf_tracer_solution(z, C::Number, ΔC::Number, profile_function::Erf)
Solution to the heat equation for a tracer concentration field `C` subject to initial
conditions that are a Heaviside (or modified Heaviside) step function aat time `t`.

## Function arguments:

- `z` for the Oceananigans model grid to evaulate the function at;
- `C` tracer value in deeper part of the step;
- `ΔC` difference in tracer between the steps;
- `κ` the diffusivity of the tracer;
- `profile_function::Erf` container with the `interface_location` and time `t` at which to
evaluate the error function solution

Default behaviour puts the `interface_location` in the centre of the depth range given by `z`.
"""
erf_tracer_solution(z, Cₗ::Number, ΔC::Number, κ::Number, profile_function::Erf) =
    Cₗ + 0.5 * ΔC * (1 + erf((z - profile_function.interface_location) /
                              sqrt(4 * κ * profile_function.time)))
"""
    tanh_initial_condition(z, Cᵤ::Number, ΔC::Number, interface_location)
Set a hyperbolic tangent initial condition for a tracer `C` over the vertical domain `z`.

## Function arguments

- `z` for the Oceananigans model grid to evaluate the function at;
- `Cˡ` tracer value in the lower layer;
- `ΔC` difference in tracer between upper layer and lower layer;
- `profile_function::HyperbolicTangent` container with `interface_location` and the
`interface_transition_width`.
"""
tanh_initial_condition(z, Cˡ::Number, ΔC::Number, profile_function::HyperbolicTangent) =
    Cˡ + 0.5 * ΔC * (1  + tanh(profile_function.interface_transition_width *
                               (z - profile_function.interface_location)))
"""
    function perturb_salinity(z, salinity_perturbation::GaussianProfile)
Perturb salinity by setting a vertical Gaussian profile in the upper layer centred at `μ`
with width `σ`.
"""
function perturb_salinity(z, salinity_perturbation::GaussianProfile)

    μ, σ, scale = salinity_perturbation.μ, salinity_perturbation.σ,
                  salinity_perturbation.scale

    if z > salinity_perturbation.interface_location
        scale * exp(-(z - μ)^2 / 2*(σ)^2) / sqrt(2*π*σ^2)
    else
        0
    end

end
"""
    function perturb_salinity(x, y, z, salinity_perturbation::GaussianBlob)
Perturb salinity by setting a horizontal Gaussian blob in the upper layer centred at `μ`
with width `σ` at depth `depth`. **Note** the `depth` needs to be an exact match to a
depth at that the `Center` in the `z` direction.
"""
function perturb_salinity(x, y, z, salinity_perturbation::GaussianBlob)

    μ, σ, scale = salinity_perturbation.μ, salinity_perturbation.σ,
                  salinity_perturbation.scale

    if z == salinity_perturbation.depth
        scale * exp(- ((x - μ[1])^2 + (y - μ[2])^2) / 2*σ^2) / (2*π*σ^2)
    else
        0
    end

end
"""
    function perturb_salinity(z, salinity_perturbation::RandomPerturbations)
Perturb salinity by adding `scale`d random noise to the salinity at `depth`.
**Note** the `depth` needs to be an exact match to a depth at that the `Center` in the `z`
direction.
"""
function perturb_salinity(z, salinity_perturbation::RandomPerturbations)

    if z == salinity_perturbation.depth
        salinity_perturbation.scale * randn()
    else
        0
    end

end
"""
    function non_dimensional_numbers(model::Oceananigans.AbstractModel,
                                     initial_conditions::TwoLayerInitialConditions)
Compute non-dimensional numbers related to the DNS experiments. The non-dimensional numbers
are:

- Prandtl number: ``Pr = ν / κₜ``
- Schmidt number: ``Sc = ν / κₛ``
- Lewis number:   ``Le = κₜ / κₛ``
- Raleigh number (density): ``Ra_{d} = Ra_{t} / Ra_{s} = (αΔT / βΔS) * (1 / Le)``.

These numbers are then saved into the simulation output file.
"""
function non_dimensional_numbers(model::Oceananigans.AbstractModel,
                                 initial_conditions::TwoLayerInitialConditions)

    ν = model.closure.ν
    κₛ, κₜ = model.closure.κ
    Pr = ν / κₜ
    Sc = ν / κₛ
    Le = κₜ / κₛ
    α = gsw_alpha(initial_conditions.S₀ˡ, initial_conditions.T₀ˡ, 0)
    β = gsw_beta(initial_conditions.S₀ˡ, initial_conditions.T₀ˡ, 0)
    Ra = ((α * initial_conditions.ΔT₀ )/ (β * initial_conditions.ΔS₀)) * (1 / Le)

    return Dict("Pr" => Pr, "Sc" => Sc, "Le" => Le, "Ra_ρ" => Ra)

end
"""
    function DNS_simulation_setup(model::Oceananigans.AbstractModel, Δt::Number,
                                  stop_time::Number,
                                  initial_conditions::TwoLayerInitialConditions)
Setup a DNS from `initial_conditions` that are of type `TwoLayerInitialConditions`.
Important non-dimensional numnbers that are part of this experiment are computed and saved
to the simulation output file.
"""
function DNCS.DNS_simulation_setup(model::Oceananigans.AbstractModel, Δt::Number,
                                   stop_time::Number, save_schedule::Number,
                                   initial_conditions::TwoLayerInitialConditions;
                                   cfl = 0.75,
                                   diffusive_cfl = 0.75,
                                   max_change = 1.2,
                                   max_Δt = 1e-1)

    simulation = Simulation(model; Δt, stop_time)

    # time step adjustments
    wizard = TimeStepWizard(; cfl, diffusive_cfl, max_change, max_Δt)
    simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

    # save output
    ϵ = KineticEnergyDissipationRate(model)
    outputs = (S = model.tracers.S, T = model.tracers.T, ϵ = ϵ, w = model.velocities.w)
    filename = form_filename(initial_conditions)
    simulation.output_writers[:outputs] = JLD2OutputWriter(model, outputs,
                                                    filename = filename,
                                                    schedule = TimeInterval(save_schedule),
                                                    overwrite_existing = true)
    jldopen(filename, "a+") do file
        file["Non_dimensional_numbers"] = non_dimensional_numbers(model, initial_conditions)
    end

    # progress reporting
    simulation.callbacks[:progress] = Callback(simulation_progress, IterationInterval(100))

    return simulation

end
"""
    function form_filename(initial_conditions::TwoLayerInitialConditions)
Create a directory based on the temperature of the upper layer and a file for the saved
output based on the type of initial condition (i.e. stable, cabbeling or unstable).
"""
function form_filename(initial_conditions::TwoLayerInitialConditions)

    ic_type = typeof(initial_conditions)
    savefile = ic_type <: StableTwoLayerInitialConditions ? "stable" :
                            ic_type <: CabbelingTwoLayerInitialConditions ?
                                "cabbeling" : ic_type <: UnstableTwoLayerInitialConditions ?
                                              "unstable" : ic_type <: IsohalineTwoLayerInitialConditions ?
                                                            "isohaline" : "isothermal"
    # make a simulation directory if one is not present
    if !isdir(SIMULATION_PATH)
        mkdir(SIMULATION_PATH)
    end
    filename = joinpath(SIMULATION_PATH, savefile * ".jld2")

    return filename

end
"""
    function compute_density(S_timeseries::FieldTimeSeries, T_timeseries::FieldTimeSeries;
                             reference_pressure = 0)
Return a density `FieldTimeSeries` calculated from the salinity and temperature
`FieldTimeSeries` DNS simulation output. The keyword argument `reference_pressure` can be
passed to specify a reference pressure at which to compute the density variable.
"""
function compute_density(S_timeseries::FieldTimeSeries, T_timeseries::FieldTimeSeries;
                         reference_pressure = 0)

    ρ_ts = FieldTimeSeries{Center, Center, Center}(S_timeseries.grid, S_timeseries.times,
                                                  indices = S_timeseries.indices,
                                                  boundary_conditions =
                                                    S_timeseries.boundary_conditions)

    t = S_timeseries.times
    for i ∈ eachindex(t)
        Sᵢ, Θᵢ = S_timeseries[i], T_timeseries[i]
        ρ_ts[i] .= @at (Center, Center, Center) gsw_rho.(Sᵢ, Θᵢ, reference_pressure)
    end

    return ρ_ts

end
"""
    function compute_density!(filepath::String; reference_pressure = 0)
Compute a density variable at `reference_pressure` and append it to saved output. **Note**
this only needs to be done once! After that the group will exist and will not be able to be
overwritten but a different group can be made by passing a different `density_string`.
This allows calling `FieldTimeSeries` on the `density_variable`
"""
function compute_density!(filepath::String; reference_pressure = 0, density_string = "ρ")

    file = jldopen(filepath, "a+")
    file_keys = keys(file["timeseries"]["S"])
    JLD2.Group(file, "timeseries/"*density_string)
    for (i, key) ∈ enumerate(file_keys)
        if i == 1
            for k ∈ keys(file["timeseries/S/"*key])
                file["timeseries/"*density_string*"/"*key*"/"*k] =
                    file["timeseries/S/"*key*"/"*k]
            end
            file["timeseries/"*density_string*"/"*key*"/reference_pressure"] =
                reference_pressure
        else
            Sᵢ, Θᵢ = file["timeseries/S/"*key], file["timeseries/T/"*key]
            file["timeseries/"*density_string*"/"*key]= gsw_rho.(Sᵢ, Θᵢ, reference_pressure)
        end
    end

    close(file)

    return nothing

end

end # module
