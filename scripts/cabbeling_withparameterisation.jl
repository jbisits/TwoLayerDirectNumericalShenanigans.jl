# Testing `quasiDNS_cabbeling`

using DirectNumericalCabbelingShenanigans

resolution = (Nx = 10, Ny = 10 , Nz = 1000)
diffusivities = (ν = 1e-6, κ = (S = 1e-9, T = 1e-7))
## Initial conditions, interface in middle of domain
S₀ = (upper = 34.568, lower = 34.7)
Θ₀ = (upper = -1.5, lower = 0.5)
reference_density = gsw_rho(S₀.lower, Θ₀.lower, 0)

model = DNS_cabbeling(resolution, diffusivities; reference_density)
model.grid

set_two_layer_initial_conditions!(model, S₀, Θ₀)

## This will be new way to set initial conditions with tanh
# Temperature
test_array_Θ = Array{Float64}(undef, size(model.grid))
Θ₀_test = range(-2, 1; length = resolution.Nz)
# ΔΘ = -2
# scale = ΔΘ / 2
# xy_centre = ΔΘ == -1 ? 0 : ΔΘ < -1 ? scale / 2 : -scale
#ΔΘ = (Θ₀.upper - Θ₀.lower) / 2
ΔΘ = (-1 - Θ₀.lower) / 2
z_centre = 0.5
interface = 10
test_initial_profile = ΔΘ .* tanh.(interface .* (Θ₀_test .+ z_centre)) .+ (Θ₀.lower + ΔΘ)
scatterlines(test_initial_profile, -1:1/resolution.Nz:0-1/resolution.Nz)

for i ∈ axes(test_array, 1), j ∈ axes(test_array, 2)
    test_array_Θ[i, j, :] = test_initial_profile
end

# Salt
test_array_S = Array{Float64}(undef, size(model.grid))
S₀_test = range(-1, 0; length = resolution.Nz)
ΔS = (S₀.upper - S₀.lower) / 2
z_centre = 0.5
interface = 10
test_initial_profile = ΔS .* tanh.(interface .* (z .+ z_centre)) .+ (S₀.lower + ΔS)
scatterlines(test_initial_profile, -1:1/resolution.Nz:0-1/resolution.Nz)

## visualise the temperature initial condition on x-z plane

x, y, z = nodes(model.grid, (Center(), Center(), Center()))
fig, ax, hm = heatmap(x, z, interior(model.tracers.S, :, 1, :); colormap = :haline)
fig, ax, plt = lines(interior(model.tracers.S, 1, 1, :), z)
Colorbar(fig[1, 2], hm)
fig

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
                                                filename = joinpath("data/simulations",
                                                                    "cabbeling.jld2"),
                                                schedule = IterationInterval(50),
                                                overwrite_existing = true)

## Progress reporting
simulation.callbacks[:progress] = Callback(simulation_progress, IterationInterval(50))

## Run the simulation
run!(simulation)
