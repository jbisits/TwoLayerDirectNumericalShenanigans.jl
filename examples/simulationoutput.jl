## Example script for viewing output
using Oceananigans.Fields
using DirectNumericalCabbelingShenanigans
using DirectNumericalCabbelingShenanigans.OutputUtilities

## Load in saved output
sim_path = joinpath(SIMULATION_PATH, "saved_output.jld2")
T_ts = FieldTimeSeries(sim_path, "T")
S_ts = FieldTimeSeries(sim_path, "S")

## Initial snapshots
visualise_snapshot(T_ts, "Θ (°C)", 1, 1, 1)
visualise_snapshot(S_ts, "S (gkg⁻ꜝ)", 1, 1, 1; colormap = :haline)

## Animations (x-z)
animate_2D_field(T_ts, "Θ (°C)", 1, 1)

## Compute a density `FieldTimeSeries`
σ₀_ts = compute_density(S_ts, T_ts)
visualise_snapshot(σ₀_ts, "σ₀ (kgm⁻³)", 1, 1, 1; colormap = :dense)
animate_2D_field(σ₀_ts, "σ₀ (kgm⁻³)", 1, 1; colorrange = (1027.68, 1027.709),
                                            colormap = cgrad(:dense)[2:end-1],
                                            lowclip = cgrad(:dense)[1],
                                            highclip = cgrad(:dense)[end])
