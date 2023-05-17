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

## Closure, diffusivity to start with parameterise a constant vertical diffusivity. Likely
#  it is better to build the simulation up from Prandtl and Raleigh numbers to give better
#  indication of what the controlling dynamics are.
ν, κ = 1e-6, 1e-5
Pr = ν / κ
closure = ScalarDiffusivity(; ν, κ)

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
        S₀[:, :, i] .= 34.7
        T₀[:, :, i] .= 0.5
    else
        S₀[:, :, i] .= 34.568
        T₀[:, :, i] .= -1.5
    end

end
set!(model, T = T₀, S = S₀)

## visualise the initial condition
x, y, z = nodes(model.grid, (Center(), Center(), Center()))
fig, ax, hm = heatmap(x, z, interior(model.tracers.T, :, 1, :))
Colorbar(fig[1, 2], hm)
fig
## simulation
Δt = 1e-2
stop_iteration = 1000
simulation = Simulation(model; Δt, stop_iteration)

## time step adjustments
wizard = TimeStepWizard(cfl = 0.75, max_Δt = 0.1, max_change = 1.2)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

## save info
outputs = (T = model.tracers.T, S = model.tracers.S)
simulation.output_writers[:outputs] = JLD2OutputWriter(model, outputs,
                                                filename = joinpath("data/simulations", "cabbeling.jld2"),
                                                schedule = IterationInterval(50))
## Progress reporting
progress(sim) = @printf("i: % 6d, sim time: % 1.3f, wall time: % 10s, Δt: % 1.4f, advective CFL: %.2e, diffusive CFL: %.2e\n",
                        iteration(sim), time(sim), prettytime(sim.run_wall_time),
                        sim.Δt, AdvectiveCFL(sim.Δt)(sim.model), DiffusiveCFL(sim.Δt)(sim.model))
simulation.callbacks[:progress] = Callback(progress, IterationInterval(50))
## Run the simulation
run!(simulation)

## Saved uutput
using Oceananigans.Fields
sim_path = joinpath(SIMULATION_PATH, "cabbeling.jld2")
Θ_ts = FieldTimeSeries(sim_path, "T")
S_ts = FieldTimeSeries(sim_path, "S")
t = Θ_ts.times

fig, ax, hm = heatmap(x, z, interior(Θ_ts, :, 1, :, 1))
Colorbar(fig[1, 2], hm)
fig

## Animations

n = Observable(1)
Θₙ = @lift interior(Θ_ts[$n], :, 1, :)
title = @lift @sprintf("t=%1.2f", t[$n])
fig, ax, hm = heatmap(x, z, Θₙ)
ax.xlabel = "x"
ax.ylabel = "z"
Colorbar(fig[1, 2], hm, label = "Temperature (°C)")
fig

frames = eachindex(t)

record(fig, "xz_temperature.mp4", frames, framerate=8) do i
    msg = string("Plotting frame ", i, " of ", frames[end])
    print(msg * " \r")
    n[] = i
end

## Density computation
σ₀_ts = deepcopy(S_ts)
for i ∈ eachindex(t)
    Sᵢ, Θᵢ = S_ts[i], Θ_ts[i]
    σ₀_ts[i] .= @at (Center, Center, Center) gsw_sigma0.(Sᵢ, Θᵢ)
end
σ₀_ts

n = Observable(1)
σ₀ⁿ = @lift interior(σ₀_ts[$n], :, 1, :)
title = @lift @sprintf("t=%1.2f", t[$n])
fig, ax, hm = heatmap(x, z, σ₀ⁿ; colormap = :dense)
ax.xlabel = "x"
ax.ylabel = "z"
Colorbar(fig[1, 2], hm, label = "σ₀ (kgm⁻³)")
fig

frames = eachindex(t)

record(fig, "xz_sigma0.mp4", frames, framerate=8) do i
    msg = string("Plotting frame ", i, " of ", frames[end])
    print(msg * " \r")
    n[] = i
end
