"""
Set up a quasi DNS (with diffusivity parameterised) as a proof-of-concept for some of the
ideas we have been exploring. The 3D domain is rectangular, periodic x-y, and no flux at
the top faces. Exploring some options which will then be moved into some execution scripts
as part of this unregistered package.
"""

using DirectNumericalShenanigans

## grid, setup is 0.5m × 0.5m × 1m currently with resolution of 1cm (might be too big),
#  variable spacing may also be an option
Lx, Ly, Lz = 0.5, 0.5, 1
Nx, Ny, Nz = 50, 50, 100
grid = RectilinearGrid(topology = (Periodic, Periodic, Bounded), size = (Nx, Ny, Nz),
                       x = (-Lx/2, Lx/2),
                       y = (-Ly/2, Ly/2),
                       z = (-Lz, 0))

## buoyancy, for now use the cabbeling EOS
eos = SeawaterPolynomials.RoquetEquationOfState(:Cabbeling)
buoyancy = SeawaterBuoyancy(equation_of_state = eos)
tracers = (:S, :T)
## By default Bounded domains have zero flux boundary conditions so leave that for now

## Closure, diffusivity to start with parameterise a constant vertical diffusivity.
ν, κ = 1e-6, 1e-7
closure = ScalarDiffusivity(; ν, κ)

## timestepper
timestepper = :RungeKutta3

## advection scheme
advection = WENO()

## model
model = NonhydrostaticModel(; grid, buoyancy, tracers, closure, timestepper, advection)

## Initial conditions, interface in middle of domain
S₀ = (upper = 34.6, lower = 34.7)
Θ₀ = (upper = -1.5, lower = 0.5)

set_two_layer_initial_conditions!(model, S₀, Θ₀)

## visualise the temperature initial condition on x-z plane
x, y, z = nodes(model.grid, (Center(), Center(), Center()))
fig, ax, hm = heatmap(x, z, interior(model.tracers.S, :, 1, :))
Colorbar(fig[1, 2], hm)
fig

## simulation
Δt = 1e-2
stop_iteration = 250
simulation = Simulation(model; Δt, stop_iteration)

## time step adjustments
wizard = TimeStepWizard(cfl = 0.75, max_Δt = 0.1, max_change = 1.2)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

## save info
outputs = (T = model.tracers.T, S = model.tracers.S)
simulation.output_writers[:outputs] = JLD2OutputWriter(model, outputs,
                                                filename = joinpath("data/simulations", "unstable.jld2"),
                                                schedule = IterationInterval(50))
## Progress reporting
progress(sim) = @printf("i: % 6d, sim time: % 1.3f, wall time: % 10s, Δt: % 1.4f, advective CFL: %.2e, diffusive CFL: %.2e\n",
                        iteration(sim), time(sim), prettytime(sim.run_wall_time),
                        sim.Δt, AdvectiveCFL(sim.Δt)(sim.model), DiffusiveCFL(sim.Δt)(sim.model))
simulation.callbacks[:progress] = Callback(progress, IterationInterval(50))
## Run the simulation
run!(simulation)
