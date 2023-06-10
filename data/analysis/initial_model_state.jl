## Note these visualisations come from a `model`, not the saved model output.

## viualise vertical grid resolution
fig, ax = scatterlines(zspacings(model.grid, Center()), znodes(model.grid, Center());
                       markersize = 7)
ax.xlabel = "Δz"
ax.ylabel = "z"
ax.title = "Vertical grid spacing for DNS"
fig

visualise_initial_conditions(model)

## visualise density profile and in x-z
σ₀ = gsw_rho.(interior(model.tracers.S, :, 1, :, 1), interior(model.tracers.T, :, 1, :), 0)
fig, ax = lines(σ₀[1, :], z)
ax.title = "Initial density profile"
ax.xlabel = "σ₀ (kgm⁻³)"
ax.ylabel = "z (m)"
fig
##
fig, ax, hm = heatmap(x, z, σ₀; colormap = :dense)
ax.title = "Initial density (x-z)"
ax.xlabel = "x (m)"
ax.ylabel = "z (m)"
Colorbar(fig[1, 2], hm)
fig
