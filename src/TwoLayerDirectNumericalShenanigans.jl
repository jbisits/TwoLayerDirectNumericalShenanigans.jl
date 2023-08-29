module TwoLayerDirectNumericalShenanigans

using Oceananigans, Printf, Reexport, JLD2, GibbsSeaWater
using SeawaterPolynomials: TEOS10EquationOfState
using Oceananigans: AbstractModel
using SpecialFunctions: erf
using Oceanostics: KineticEnergyDissipationRate
import Base: show

@reexport using Oceananigans, Reexport

export TwoLayerDNS, DNS, DNS_simulation_setup

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
    IsothermalTwoLayerInitialConditions

export set_two_layer_initial_conditions!, add_velocity_random_noise!

export ContinuousProfileFunction, HyperbolicTangent, Erf

export SalinityPerturbation, GaussianProfile, GaussianBlob, RandomPerturbations

export DOMAIN_EXTENT, HIGH_RESOLUTION, SO_DIFFUSIVITIES, REFERENCE_DENSITY,
       INTERFACE_LOCATION, SIMULATION_PATH, DNCS

export compute_density, compute_density!, non_dimensional_numbers

export animate_2D_field, visualise_initial_conditions, visualise_initial_density,
       visualise_snapshot

include("initialconditions.jl")
include("continuousprofilefunctions.jl")
include("salinityperturbations.jl")
include("twolayerdns.jl")
include("set_initialconditions.jl")
include("makiefunctions.jl")
include("constants.jl")
include("output_utils.jl")

end
