# Testing `quasiDNS_cabbeling`

using DirectNumericalShenanigans

resolution = (Nx = 10, Ny = 10 , Nz = 1000)
diffusivities = (ν = 1e-6, κ = (S = 1e-9, T = 1e-7))
## Initial conditions, interface in middle of domain
S₀ = (upper = 34.568, lower = 34.7)
Θ₀ = (upper = -1.5, lower = 0.5)
reference_density = gsw_rho(S₀.lower, Θ₀.lower, 0)

model = quasiDNS_cabbeling(resolution, diffusivities; reference_density)
model.grid
set_two_layer_initial_conditions!(model, S₀, Θ₀)

## This will be new way to set initial conditions with tanh
test_array = size(model.grid)
Θ₀_test = range(Θ₀.upper, Θ₀.lower; length = resolution.Nz)

scatterlines(-tanh.(10 .* (Θ₀_test .+ 0.5)), -1:1/resolution.Nz:0-1/resolution.Nz)

## visualise the temperature initial condition on x-z plane
x, y, z = nodes(model.grid, (Center(), Center(), Center()))
fig, ax, hm = heatmap(x, z, interior(model.tracers.S, :, 1, :); colormap = :haline)
# Colorbar(fig[1, 2], hm)
# fig

## Random noise in horizontal velocities
u, v, w = model.velocities
ϵ = 1e-4
uᵢ, vᵢ = ϵ * randn(size(u)), ϵ * randn(size(v))
set!(model, u = uᵢ, v = vᵢ)

## simulation
Δt = 1e-2
stop_iteration = 100
simulation = Simulation(model; Δt, stop_iteration)

## time step adjustments
wizard = TimeStepWizard(cfl = 0.75, max_Δt = 0.1, max_change = 1.2)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

## save info
outputs = (T = model.tracers.T, S = model.tracers.S)
simulation.output_writers[:outputs] = JLD2OutputWriter(model, outputs,
                                                filename = joinpath("data/simulations", "cabbeling.jld2"),
                                                schedule = IterationInterval(50),
                                                overwrite_existing = true)

## Progress reporting
simulation.callbacks[:progress] = Callback(simulation_progress, IterationInterval(50))

## Run the simulation
run!(simulation)
