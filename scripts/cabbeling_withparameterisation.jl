"""
Set up a quasi DNS (with diffusivity parameterised) as a proof-of-concept for some of the
ideas we have been exploring. The 3D domain is rectangular, periodic x-y, and no flux at
the top faces. Exploring some options which will then be moved into some execution scripts
as part of this unregistered package.
"""

using DirectNumericalShenanigans

## grid, setup is 0.5m × 0.5m × 1m currently with resolution of 10cm (might be too big),
#  variable spacing may also be an option
x_length, y_length, z_length = 0.5, 0.5, 1
Nx, Ny, Nz = 5, 5, 10
grid = RectilinearGrid(topology = (Periodic, Periodic, Bounded),
                       extent   = (x_length, y_length, z_length),
                       size     = (Nx, Ny, Nz))

## buoyancy, for now use the cabbeling EOS
eos = SeawaterPolynomials.RoquetEquationOfState(:Cabbeling)
buoyancy = SeawaterBuoyancy(equation_of_state = eos)
tracers = (:S, :T)
## By default Bounded domains have zero flux boundary conditions so leave that for now

## Closure, diffusivity to start with parameterise a constant vertical diffusivity
ν, κ = 1e-6, 1e-5
closure = ScalarDiffusivity(; ν, κ)

## model
model = NonhydrostaticModel(; grid, buoyancy, tracers, closure, timestepper = :RungeKutta3)

## Initial conditions, interface in middle of domain
S₀ = Array{Float64}(undef, size(model.grid))
T₀ = Array{Float64}(undef, size(model.grid))

for i ∈ 1:10

    if i ≤ floor(Nz) / 2
        S₀[:, :, i] .= 34.568
        T₀[:, :, i] .= -1.5
    else
        S₀[:, :, i] .= 34.7
        T₀[:, :, i] .= 0.5
    end

end

set!(model, T = T₀, S = S₀)

## visualise the initial condition
interior(model.tracers.T, :, :, :)
heatmap(interior(model.tracers.T, :, 1, :))

## simulation
Δt = κ / (z_length / Nz)
stop_iteration = 100
simulation = Simulation(model; Δt, stop_iteration)

## save info
outputs = (T = model.tracers.T, S = model.tracers.S)
simulation.output_writers[:outputs] = JLD2OutputWriter(model, outputs,
                                                filename = joinpath("data/simulations", "cabbeling.jld2"),
                                                schedule = TimeInterval(Δt * 10))
