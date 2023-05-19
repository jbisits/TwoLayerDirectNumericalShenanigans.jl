## Saved Output
using Oceananigans.Fields, GibbsSeaWater

sim_path = joinpath(SIMULATION_PATH, "unstable.jld2")
Θ_ts = FieldTimeSeries(sim_path, "T")
S_ts = FieldTimeSeries(sim_path, "S")
x, y, z = nodes(model.grid, (Center(), Center(), Center()))
t = Θ_ts.times

σ₀_ts = deepcopy(S_ts)
for i ∈ eachindex(t)
    Sᵢ, Θᵢ = S_ts[i], Θ_ts[i]
    σ₀_ts[i] .= @at (Center, Center, Center) gsw_sigma0.(Sᵢ, Θᵢ)
end
σ₀_ts

## Plots (x-z)
fig, ax, hm = heatmap(x, z, interior(Θ_ts, :, 1, :, 1); colormap = :thermal)
Colorbar(fig[1, 2], hm)
fig

fig, ax, hm = heatmap(x, z, interior(σ₀_ts, :, 1, :, 1); colormap = :dense)
Colorbar(fig[1, 2], hm)
fig

# Animations (x-z)

## Temperature
n = Observable(1)
Θₙ = @lift interior(Θ_ts[$n], :, 1, :)
title = @lift @sprintf("t=%1.2f", t[$n])
fig, ax, hm = heatmap(x, z, Θₙ; colormap = :thermal)
ax.xlabel = "x"
ax.ylabel = "z"
Colorbar(fig[1, 2], hm, label = "Temperature (°C)")
fig

frames = eachindex(t)

record(fig, "xz_temperature_unstable.mp4", frames, framerate=8) do i
    msg = string("Plotting frame ", i, " of ", frames[end])
    print(msg * " \r")
    n[] = i
end

## Density
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
