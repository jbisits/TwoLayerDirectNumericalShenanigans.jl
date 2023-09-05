# A two layer DNS experiment example
using TwoLayerDirectNumericalShenanigans

architecture = CPU() # or GPU()
diffusivities = (ν = 1e-4, κ = (S = 1e-5, T = 1e-5))
resolution = (Nx = 10, Ny = 10, Nz = 100)

## Setup the model
model = DNS(architecture, DOMAIN_EXTENT, resolution, diffusivities;
            reference_density = REFERENCE_DENSITY)

## set initial conditions, currently there are four options available in this submodule
T₀ᵘ = -1.5
S₀ᵘ = (stable = 34.551, cabbeling = 34.568, unstable = 34.59)
stable = StableUpperLayerInitialConditions(S₀ᵘ.stable, T₀ᵘ)
initial_conditions = TwoLayerInitialConditions(stable)
profile_function = HyperbolicTangent(INTERFACE_LOCATION, 50.0)
tracer_perturbation_depth = find_depth(model, INTERFACE_LOCATION / 2)
tracer_perturbation = SalinityGaussianProfile(tracer_perturbation_depth, 0.0, 1.5)
noise_depth = find_depth(model, INTERFACE_LOCATION)
initial_noise = SalinityNoise(noise_depth, 2.0)

dns = TwoLayerDNS(model, profile_function, initial_conditions; tracer_perturbation, initial_noise)

set_two_layer_initial_conditions!(dns)
# using CairoMakie - need to have a different environment activated with CairoMakie as dep
# visualise_initial_conditions(dns, 1, 1)
# visualise_initial_density(dns, 1, 1, 0)

## build the simulation
Δt = 1e-4
stop_time = 60
save_schedule = 0.5 # seconds
simulation = DNS_simulation_setup(dns, Δt, stop_time, save_schedule)

## Run the simulation
run!(simulation)
