## File to view saved output
using Oceananigans.Fields
using DirectNumericalCabbelingShenanigans
using DirectNumericalCabbelingShenanigans.OutputUtilities

## Load in saved output
sim_path = joinpath(SIMULATION_PATH, "cabbeling.jld2")
T_ts = FieldTimeSeries(sim_path, "T")
S_ts = FieldTimeSeries(sim_path, "S")

## Snapshots
visualise_snapshot(T_ts, "Θ (°C)", 6)
visualise_snapshot(S_ts, "S (gkg⁻ꜝ)", 2; colormap = :haline)

## Animations (x-z)
animate_2D_field(T_ts, "Θ (°C)", (x = x, z = z))

## Compute a density `FieldTimeSeries`
σ₀_ts = compute_density(S_ts, T_ts)
visualise_snapshot(σ₀_ts, "σ₀ (kgm⁻³)", 6; colormap = :dense)
animate_2D_field(σ₀_ts, "σ₀ (kgm⁻³)", (x = x, z = z))
