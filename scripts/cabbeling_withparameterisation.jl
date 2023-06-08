# Testing `quasiDNS_cabbeling`

using DirectNumericalCabbelingShenanigans

architecture = CPU() # or GPU()
resolution = (Nx = 50, Ny = 50, Nz = 1000)
diffusivities = (ν = 1e-4, κ = (S = 1e-5, T = 1e-5))
## Initial conditions
S₀ = (upper = 34.6, lower = 34.7)
Θ₀ = (upper = -1.5, lower = 0.5)
reference_density = gsw_rho(S₀.lower, Θ₀.lower, 0)

model = DNS_cabbeling(architecture, resolution, diffusivities; reference_density)
model.grid

## viualise vertical grid resolution
fig, ax = scatterlines(zspacings(model.grid, Center()), znodes(model.grid, Center());
                       markersize = 7)
ax.xlabel = "Δz"
ax.ylabel = "z"
ax.title = "Vertical grid spacing for DNS"
fig

## set initial conditions
set_two_layer_initial_conditions!(model, S₀, Θ₀;
                                  interface_location = 0.375, interface_thickness = 100,
                                  salinity_perturbation_width = 200)
non_dimensional_numbers(model, S₀, Θ₀)
## visualise the salt initial condition on x-z plane
x, y, z = nodes(model.grid, (Center(), Center(), Center()))
fig, ax, hm = heatmap(x, z, interior(model.tracers.S, :, 1, :); colormap = :haline)
ax.title = "Initial salinity (x-z)"
ax.xlabel = "x (m)"
ax.ylabel = "z (m)"
Colorbar(fig[1, 2], hm)
fig

## visualise the salinity profile
fig, ax = lines(interior(model.tracers.S, 1, 1, :), z)
ax.title = "Initial salinity profile"
ax.xlabel = "S (gkg⁻¹)"
ax.ylabel = "z (m)"
fig
## visualise temperature profile
fig, ax = lines(interior(model.tracers.T, 1, 1, :), z)
ax.title = "Initial temperature profile"
ax.xlabel = "Θ (°C)"
ax.ylabel = "z (m)"
fig
## visualise density profile and in x-z
σ₀ = gsw_rho.(interior(model.tracers.S, :, 1, :), interior(model.tracers.T, :, 1, :), 0)
fig, ax = lines(σ₀[1, :], z)
ax.title = "Initial density profile"
ax.xlabel = "σ₀ (kgm⁻³)"
ax.ylabel = "z (m)"
fig
##
fig, ax, hm = heatmap(x, z, σ₀; colormap = :dense)
ax.title = "Initial density (x-z)"
ax.xlabel = "x (m)"
ax.ylabel = "z (m)"
Colorbar(fig[1, 2], hm)
fig

## simulation
Δt = 1e-4
stop_iteration = 100
simulation = Simulation(model; Δt, stop_iteration)

## time step adjustments
wizard = TimeStepWizard(cfl = 0.75, max_Δt = 1e-2, max_change = 1.2)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

## save info
outputs = (S = model.tracers.S, T = model.tracers.T)
simulation.output_writers[:outputs] = JLD2OutputWriter(model, outputs,
                                                filename = joinpath("data/simulations",
                                                                    "unstable.jld2"),
                                                schedule = IterationInterval(50),
                                                overwrite_existing = true)

## Progress reporting
simulation.callbacks[:progress] = Callback(simulation_progress, IterationInterval(50))

## Run the simulation
run!(simulation)
