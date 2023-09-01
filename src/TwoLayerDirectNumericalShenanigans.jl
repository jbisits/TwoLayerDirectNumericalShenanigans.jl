module TwoLayerDirectNumericalShenanigans

using Oceananigans, Printf, Reexport, JLD2, GibbsSeaWater
using SeawaterPolynomials: TEOS10EquationOfState
using Oceananigans: AbstractModel
using SpecialFunctions: erf
using Oceanostics: KineticEnergyDissipationRate
using CUDA: allowscalar
import Base: show, iterate

@reexport using Oceananigans, Reexport

"Abstract type for `TwoLayerDNS`."
abstract type AbstractTwoLayerDNS end
"Abstract super type for a perturbation added to a tracer (`S` or `T`)."
abstract type AbstractTracerPerturbation end
"Abstract super type for random noise added to tracer or velocity field."
abstract type AbstractNoise end
"Abstract super type for the continuous function that sets the continuous profile for
temperature and salinity."
abstract type AbstractContinuousProfileFunction end
"Abstract supertype for dns initial conditions."
abstract type AbstractInitialConditions end

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

export set_two_layer_initial_conditions!, find_depth

export AbstractContinuousProfileFunction, HyperbolicTangent, Erf

export AbstractTracerPerturbation, SalinityGaussianProfile, SalinityGaussianBlob,
       TemperatureGaussianProfile, TemperatureGaussianBlob

export AbstractNoise, SalinityNoise, TemperatureNoise, VelocityNoise

export DOMAIN_EXTENT, HIGH_RESOLUTION, SO_DIFFUSIVITIES, REFERENCE_DENSITY,
       INTERFACE_LOCATION, SIMULATION_PATH, TLDNS

export compute_density, compute_density!, non_dimensional_numbers

export animate_2D_field, visualise_initial_conditions, visualise_initial_density,
       visualise_snapshot

include("initialconditions.jl")
include("continuousprofilefunctions.jl")
include("tracerperturbations.jl")
include("twolayerdns.jl")
include("noiseperturbations.jl")
include("set_initialconditions.jl")
include("makiefunctions.jl")
include("constants.jl")
include("output_utils.jl")

end
