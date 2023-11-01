@inline tracer_perturbationᶜᶜᶜ(i, j, k, grid, tracer, tracer_mean) =
    tracer[i, j, k] - tracer_mean[i, j, k]
# Only needed for testing
@inline function tracer_perturbation(model, tracer)

    grid = model.grid
    tracer_mean = Field(Average(tracer))

    return KernelFunctionOperation{Center, Center, Center}(tracer_perturbationᶜᶜᶜ, grid, tracer, tracer_mean)
end

@inline vtf(i, j, k, grid, tracer, tracer_mean, w) =
        -ℑzᵃᵃᶜ(i, j, k, w) * tracer_perturbationᶜᶜᶜ(i, j, k, grid, tracer, tracer_mean)
@inline function vertical_tracer_flux(model, tracer)

    grid = model.grid
    tracer_mean = Field(Average(tracer))
    w = model.velocities.w

    return KernelFunctionOperation{Center, Center, Center}(vtf, grid, tracer, tracer_mean, w)
end
