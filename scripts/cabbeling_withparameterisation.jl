"""
Set up a quasi DNS (with diffusivity parameterised) as a proof-of-concept for some of the
ideas we have been exploring. The 3D domain is rectangular, periodic x-y, and no flux at
the top faces. Exploring some options which will then be moved into some execution scripts
as part of this unregistered package.
"""

using DirectNumericalShenanigans

## grid, setup is 0.5m × 0.5m × 1m currently with resolution of 10cm (might be too big),
#  variable spacing may also be an option
Lx, Ly, Lz = 0.5, 0.5, 1
Nx, Ny, Nz = 50, 50, 100
grid = RectilinearGrid(topology = (Periodic, Periodic, Bounded),
                       size     = (Nx, Ny, Nz),
                       x = (-Lx/2, Lx/2),
                       y = (-Ly/2, Ly/2),
                       z = (-Lz, 0))

## buoyancy, for now use the cabbeling EOS
eos = SeawaterPolynomials.RoquetEquationOfState(:Cabbeling)
buoyancy = SeawaterBuoyancy(equation_of_state = eos)
tracers = (:S, :T)
## By default Bounded domains have zero flux boundary conditions so leave that for now

## Closure, diffusivity to start with parameterise a constant vertical diffusivity
ν, κ = 1e-6, 1e-5
Pr = ν / κ
closure = VerticalScalarDiffusivity(; ν, κ)

## timestepper
timestepper = :RungeKutta3

## advection scheme
advection = WENO()

## model
model = NonhydrostaticModel(; grid, buoyancy, tracers, closure, timestepper, advection)

## Initial conditions, interface in middle of domain
S₀ = Array{Float64}(undef, size(model.grid))
T₀ = Array{Float64}(undef, size(model.grid))

for i ∈ axes(S₀, 3)

    if i ≤ floor(Nz) / 2
        S₀[:, :, i] .= 34.568
        T₀[:, :, i] .= -1.5
    else
        S₀[:, :, i] .= 34.7
        T₀[:, :, i] .= 0.5
    end

end
S₀
set!(model, T = T₀, S = S₀)

## visualise the initial condition
x, y, z = nodes(model.grid, (Center(), Center(), Center()))
heatmap(x, z, interior(model.tracers.T, :, 1, :))

## simulation
Δt = κ / (Lz / Nz)
wizard = TimeStepWizard(cfl = 0.75, max_change = 1.2; max_Δt = Δt * 10)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))
stop_iteration = 1000
simulation = Simulation(model; Δt, stop_iteration)

## save info
outputs = (T = model.tracers.T, S = model.tracers.S)
simulation.output_writers[:outputs] = JLD2OutputWriter(model, outputs,
                                                filename = joinpath("data/simulations", "cabbeling.jld2"),
                                                schedule = TimeInterval(Δt * 10))
## Progress reporting
progress(sim) = @printf("i: % 6d, sim time: % 1.3f, wall time: % 10s, Δt: % 1.4f, advective CFL: %.2e, diffusive CFL: %.2e\n",
                        iteration(sim), time(sim), prettytime(sim.run_wall_time),
                        sim.Δt, AdvectiveCFL(sim.Δt)(sim.model), DiffusiveCFL(sim.Δt)(sim.model))
simulation.callbacks[:progress] = Callback(progress, IterationInterval(50))
## Run the simulation
run!(simulation)

## Output
using Oceananigans.Fields
sim_path = joinpath(SIMULATION_PATH, "cabbeling.jld2")
Θ_ts = FieldTimeSeries(sim_path, "T")
t = Θ_ts.times
