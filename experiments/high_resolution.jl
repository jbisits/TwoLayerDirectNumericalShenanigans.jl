# High resolution model

using DirectNumericalCabbelingShenanigans

architecture = CPU() # or GPU()
diffusivities = (ν = 1e-6, κ = (S = 1e-9, T = 1e-7))
resolution = (Nx = 50, Ny = 50, Nz = 7000)

## Initial salinity and temperature conditions
upper_salinity = (stable = 34.551, cabbeling = 34.568, unstable = 34.59)
S₀ = (upper = upper_salinity.cabbeling, lower = 34.7)
Θ₀ = (upper = -1.5, lower = 0.5)
reference_density = gsw_rho(S₀.lower, Θ₀.lower, 0)

model = DNS_cabbeling(architecture, resolution, diffusivities; reference_density)

## set initial conditions
set_two_layer_initial_conditions!(model, S₀, Θ₀;
                                  interface_location = 0.375, interface_thickness = 100,
                                  salinity_perturbation_width = 100)

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
saved_ouput = joinpath(SIMULAITON_PATH, filename * ".jld2")
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
