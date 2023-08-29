# A two layer DNS experiment example
using TwoLayerDirectNumericalShenanigans

architecture = CPU() # or GPU()
diffusivities = (ν = 1e-4, κ = (S = 1e-5, T = 1e-5))

## Setup the model
model = DNS(architecture, DOMAIN_EXTENT, HIGH_RESOLUTION, diffusivities;
            reference_density = REFERENCE_DENSITY)

## set initial conditions, currently there are four options available in this submodule
T₀ᵘ = -1.5
S₀ᵘ = (stable = 34.551, cabbeling = 34.568, unstable = 34.59)
stable = StableUpperLayerInitialConditions(S₀ᵘ.stable, T₀ᵘ)
initial_conditions = TwoLayerInitialConditions(stable)
profile_function = HyperbolicTangent(INTERFACE_LOCATION, 50.0)
z = znodes(model.grid, Center(), Center(), Center())
depth_idx = findfirst(z .> 2 * INTERFACE_LOCATION / 3)
salinity_perturbation = SalinityGaussianProfile(z[depth_idx], 0.0, 1.5)
initial_noise = SalinityNoise(z[depth_idx], 2.0)

dns = TwoLayerDNS(model, profile_function, initial_conditions; salinity_perturbation, initial_noise)

set_two_layer_initial_conditions!(dns)
# using CairoMakie - need to have a different environment activated with CairoMakie as dep
# visualise_initial_conditions(dns, 1, 1)
# visualise_initial_density(dns, 1, 1, 0)

## build the simulation
Δt = 1e-4
stop_time = 1
save_schedule = 0.5 # seconds
simulation = DNS_simulation_setup(dns, Δt, stop_time, save_schedule)

## Run the simulation
run!(simulation)
