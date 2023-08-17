# Medium resolution two layer simulation
using DirectNumericalCabbelingShenanigans
using DirectNumericalCabbelingShenanigans.TwoLayerDNS

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
salinity_perturbation = GaussianBlob(z[depth_idx], [0.0, 0.0], 1.5)
set_two_layer_initial_conditions!(model, initial_conditions, profile_function,
                                  salinity_perturbation)

DNCS.OutputUtilities.visualise_initial_conditions(model, 1, 1)
DNCS.OutputUtilities.visualise_initial_density(model, 1, 1, 0)

## build the simulation
Δt = 1e-4
stop_time = 1
save_schedule = 0.5 # seconds
simulation = DNS_simulation_setup(model, Δt, stop_time, save_schedule, initial_conditions)

## Run the simulation
run!(simulation)
