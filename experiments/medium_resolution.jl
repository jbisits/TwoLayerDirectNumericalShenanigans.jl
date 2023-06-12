# Medium resolution model

using DirectNumericalCabbelingShenanigans
using DirectNumericalCabbelingShenanigans.TwoLayerDNS

architecture = CPU() # or GPU()
diffusivities = (ν = 1e-6, κ = (S = 1e-9, T = 1e-7))
domain_extent = (Lx = 0.1, Ly = 0.1, Lz = 1)
resolution = (Nx = 30, Ny = 30, Nz = 1000)
reference_density = gsw_rho(S₀ˡ, T₀ˡ, 0)

## Setup the model
model = DNS(architecture, domain_extent, resolution, diffusivities; reference_density)

## set initial conditions
T₀ᵘ = -1.5
S₀ᵘ = (stable = 34.551, cabbeling = 34.568, unstable = 34.59)
initial_conditions = CabbelingTwoLayerInitialConditions(S₀ᵘ.cabbeling, T₀ᵘ)
set_two_layer_initial_conditions!(model, initial_conditions;
                                  interface_location = 0.375, interface_thickness = 100,
                                  salinity_perturbation_width = 100)
visualise_initial_conditions(model)
## build the simulation
Δt = 1e-4
stop_time = 1
simulation = DNS_simulation_setup(model, Δt, stop_time, initial_conditions)

## Run the simulation
run!(simulation)

##
Δt = 1e-4
stop_iteration = 100
simulation = Simulation(model; Δt, stop_iteration)

## time step adjustments
wizard = TimeStepWizard(cfl = 0.75, max_Δt = 1e-2, max_change = 1.2)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

## save info
outputs = (S = model.tracers.S, T = model.tracers.T)
filename =  S₀.upper == upper_salinity.stable ? "stable" :
                        S₀.upper == upper_salinity.cabbeling ? "cabbeling" : "unstable"
saved_ouput = joinpath(SIMULATION_PATH, filename * ".jld2")
simulation.output_writers[:outputs] = JLD2OutputWriter(model, outputs,
                                                filename = saved_ouput,
                                                schedule = IterationInterval(50),
                                                overwrite_existing = true)
jldopen(saved_ouput, "a+") do file
    file["Non_dimensional_numbers"] = non_dimensional_numbers(model, S₀, Θ₀)
end

## Progress reporting
simulation.callbacks[:progress] = Callback(simulation_progress, IterationInterval(50))

## Run the simulation
run!(simulation)
