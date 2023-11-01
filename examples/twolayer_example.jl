# A two layer DNS experiment example
using TwoLayerDirectNumericalShenanigans

architecture = CPU() # or GPU()
diffusivities = (ν = 1e-4, κ = (S = 1e-5, T = 1e-5))
resolution = (Nx = 10, Ny = 10, Nz = 100)

## Setup the model
model = DNSModel(architecture, DOMAIN_EXTENT, resolution, diffusivities;
                 reference_density = REFERENCE_DENSITY)

## set initial conditions
T₀ᵘ = -1.5
S₀ᵘ = (stable = 34.551, cabbeling = 34.58, unstable = 34.59)
cabbeling = CabbelingUpperLayerInitialConditions(S₀ᵘ.cabbeling, T₀ᵘ)
initial_conditions = TwoLayerInitialConditions(cabbeling)
transition_depth = find_depth(model, INTERFACE_LOCATION)
profile_function = StepChange(transition_depth)
tracer_perturbation_depth = find_depth(model, INTERFACE_LOCATION / 2)
tracer_perturbation = SalinityGaussianProfile(tracer_perturbation_depth, 0.0, 1.5)
noise_depth = find_depth(model, INTERFACE_LOCATION)
initial_noise = SalinityNoise(noise_depth, 1e-2)

tldns = TwoLayerDNS(model, profile_function, initial_conditions; initial_noise)

set_two_layer_initial_conditions!(tldns)

## build the simulation
Δt = 1e-4
stop_time = 60
save_schedule = 0.5 # seconds
output_path = joinpath(@__DIR__, "outputs/")
checkpointer_time_interval = 30
simulation = TLDNS_simulation_setup(tldns, Δt, stop_time, save_schedule, TLDNS.save_computed_output!;
                                    output_path, checkpointer_time_interval)

## Run the simulation
run!(simulation)
