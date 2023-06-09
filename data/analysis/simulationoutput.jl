## Saved Output
using Oceananigans.Fields, GibbsSeaWater

sim_path = joinpath(SIMULATION_PATH, "cabbeling.jld2")
Θ_ts = FieldTimeSeries(sim_path, "T")
S_ts = FieldTimeSeries(sim_path, "S")
x, y, z = nodes(Θ_ts[1])
t = Θ_ts.times

σ₀_ts = deepcopy(S_ts)
for i ∈ eachindex(t)
    Sᵢ, Θᵢ = S_ts[i], Θ_ts[i]
    σ₀_ts[i] .= @at (Center, Center, Center) gsw_sigma0.(Sᵢ, Θᵢ)
end

## Plots (x-z)
fig, ax, hm = heatmap(x, z, interior(S_ts, :, 1, :, 1); colormap = :haline)
Colorbar(fig[1, 2], hm)
fig

fig, ax, hm = heatmap(x, z, interior(Θ_ts, :, 1, :, 1); colormap = :thermal)
Colorbar(fig[1, 2], hm)
fig

## Animations (x-z)
animate_2D_field(Θ_ts, "Θ (°C)", (x = x, z = z))
animate_2D_field(σ₀_ts, "σ₀ (kgm⁻³)", (x = x, z = z))
