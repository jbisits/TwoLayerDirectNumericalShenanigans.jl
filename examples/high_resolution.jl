# High resolution two layer simulation
using DirectNumericalCabbelingShenanigans
using DirectNumericalCabbelingShenanigans.TwoLayerDNS

architecture = CPU() # or GPU()
diffusivities = (ν = 1e-4, κ = (S = 1e-6, T = 1e-5))

## Setup the model
model = DNS(architecture, DOMAIN_EXTENT, HIGH_RESOLUTION, diffusivities;
            reference_density = REFERENCE_DENSITY)

## set initial conditions, currently there are four options available in this submodule
T₀ᵘ = -1.5
S₀ᵘ = (stable = 34.551, cabbeling = 34.568, unstable = 34.6)
stable = StableUpperLayerInitialConditions(S₀ᵘ.stable, T₀ᵘ)
cabbeling = CabbelingUpperLayerInitialConditions(S₀ᵘ.cabbeling, T₀ᵘ)
unstable = UnstableUpperLayerInitialConditions(S₀ᵘ.unstable, T₀ᵘ)
isohaline = IsohalineUpperLayerInitialConditions(T₀ᵘ)
initial_conditions = TwoLayerInitialConditions(unstable)
profile_function = HyperbolicTangent(INTERFACE_LOCATION, 50)
set_two_layer_initial_conditions!(model, initial_conditions, profile_function;
                                  salinity_perturbation = true)
add_velocity_random_noise!(model, 1e-2, INTERFACE_LOCATION)

## build the simulation
Δt = 1e-4
stop_time = 2
save_schedule = 0.5 # seconds
simulation = DNS_simulation_setup(model, Δt, stop_time, save_schedule, initial_conditions)

## Run the simulation
run!(simulation)
