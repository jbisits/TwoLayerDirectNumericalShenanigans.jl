@inline tracer_perturbationᶜᶜᶜ(i, j, k, grid, tracer, tracer_mean) =
    tracer[i, j, k] - tracer_mean[i, j, k]
# Only needed for testing
"""
    function tracer_perturbation(model, tracer)
Calculate the perturbation of `tracer` from the spatial `tracer` mean.
"""
@inline function tracer_perturbation(model, tracer)

    grid = model.grid
    tracer_mean = Field(Average(tracer))

    return KernelFunctionOperation{Center, Center, Center}(tracer_perturbationᶜᶜᶜ, grid, tracer, tracer_mean)
end

@inline vtf(i, j, k, grid, tracer, tracer_mean, w) =
        -ℑzᵃᵃᶜ(i, j, k, w) * tracer_perturbationᶜᶜᶜ(i, j, k, grid, tracer, tracer_mean)
"""
    function vertical_tracer_flux(model, tracer)
Calculate the vertical tracer flux of `tracer` from `model`.
"""
@inline function vertical_tracer_flux(model, tracer)

    grid = model.grid
    tracer_mean = Field(Average(tracer))
    w = model.velocities.w

    return KernelFunctionOperation{Center, Center, Center}(vtf, grid, tracer, tracer_mean, w)
end
"""
    function potential_energy(model)
Calculate the potential energy
```math
Eₚ = g∫ᵥρzdV
```
for `model` during a simulation.
"""
@inline function potential_energy(model)

    grid = model.grid
    ρ = seawater_density(model)
    Z = Oceananigans.Models.model_geopotential_height(model)
    g = tuple(model.buoyancy.model.gravitational_acceleration)

    return KernelFunctionOperation{Center, Center, Center}(Ep, grid, ρ, Z, g)
end
@inline Ep(i, j, k, grid, ρ, Z, g) = g[1] * ρ[i, j, k] * Z[i, j, k]
