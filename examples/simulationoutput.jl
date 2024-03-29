## Example script for viewing output
using TwoLayerDirectNumericalShenanigans
using CairoMakie

sim_path_nc = joinpath(SIMULATION_PATH, "stable_tanh_salinitygaussianprofile_1min.nc")
sim_path_jld2 = joinpath(SIMULATION_PATH, "stable_tanh_salinitygaussianprofile_1min.jld2")

## Load in saved output if `.nc`
using NCDatasets
ds = NCDataset(sim_path_nc)
# extract variables
close(ds)

# or using Rasters.jl
using Rasters
S_rs = Raster(sim_path_nc, name = :S)

## Load in saved output if `.jld2`
using Oceananigans.Fields
T_ts = FieldTimeSeries(sim_path, "T", backend = OnDisk())
S_ts = FieldTimeSeries(sim_path, "S", backend = OnDisk())

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
