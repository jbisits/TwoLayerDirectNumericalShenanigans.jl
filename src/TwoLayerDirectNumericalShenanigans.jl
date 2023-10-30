module TwoLayerDirectNumericalShenanigans

using Oceananigans, Printf, Reexport, JLD2, Rasters, NCDatasets, GibbsSeaWater
using Oceananigans: AbstractModel, Operators.ℑzᵃᵃᶜ
using Oceananigans: BuoyancyModels.get_temperature_and_salinity, BuoyancyModels.θ_and_sᴬ,  BuoyancyModels.Zᶜᶜᶜ
using Oceananigans: BuoyancyModels.buoyancy_perturbationᶜᶜᶜ, BuoyancyModels.∂z_b
using Oceananigans: Models.seawater_density
using SeawaterPolynomials
using SeawaterPolynomials: TEOS10EquationOfState
import SeawaterPolynomials.ρ
using SpecialFunctions: erf
using Oceanostics: KineticEnergyDissipationRate
using OceanRasterConversions: get_σₚ
using CUDA: allowscalar, CuArray
import Base: show, iterate

@reexport using Oceananigans, Reexport

"Abstract type for `TwoLayerDNS`."
abstract type AbstractTwoLayerDNS end
"Abstract super type for a perturbation added to a tracer (`S` or `T`)."
abstract type AbstractTracerPerturbation end
"Abstract super type for random noise added to tracer or velocity field."
abstract type AbstractNoise end
"Abstract super type for the function that sets the profile for
temperature and salinity."
abstract type AbstractProfileFunction end
"Abstract supertype for dns initial conditions."
abstract type AbstractInitialConditions end

export TwoLayerDNS, DNSModel, TLDNS_simulation_setup

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

export AbstractContinuousProfileFunction, HyperbolicTangent, Erf, MidPoint

export AbststractStepChangeProfileFuncion, StepChange

export AbstractTracerPerturbation, SalinityGaussianProfile, SalinityGaussianBlob,
       TemperatureGaussianProfile, TemperatureGaussianBlob

export AbstractNoise, SalinityNoise, TemperatureNoise, VelocityNoise

export DOMAIN_EXTENT, HIGH_RESOLUTION, SO_DIFFUSIVITIES, REFERENCE_DENSITY,
       INTERFACE_LOCATION, SIMULATION_PATH, CHECKPOINT_PATH, TLDNS

export compute_density, compute_density!

export animate_2D_field, visualise_initial_conditions, visualise_initial_density,
       visualise_snapshot, plot_scalar_diagnostics, hovmoller

include("initialconditions.jl")
include("continuousprofilefunctions.jl")
include("stepchangeprofilefunction.jl")
include("tracerperturbations.jl")
include("kernelfunctions.jl")
include("twolayerdns.jl")
include("noiseperturbations.jl")
include("set_initialconditions.jl")
include("makiefunctions.jl")
include("constants.jl")
include("output_utils.jl")

end
