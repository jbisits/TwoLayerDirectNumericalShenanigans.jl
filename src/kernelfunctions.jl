"""
    C(i, j, k, grid, C)
Get tracer `C` values for use in other function. There may be another way to do this for
`KernelFunctionOperation`s but I have not found it so will use this for now. **Note:** this
return the value of the tracer with no interpolation so if the tracer `C` is at
`(Center, Center, Center)` the value extracted will be at (Center, Center, Center)`.
"""
C(i, j, k, grid, C) = C[i, j, k]
"""
    σ(i, j, k, grid, model, reference_pressure)
Compute potential density `σ` at `reference_pressure` from salinity and temperature tracers
in `model`.
"""
σ(i, j, k, grid, model, reference_pressure) = gsw_rho(C(i, j, k, grid, model.tracers.S),
                                                      C(i, j, k, grid, model.tracers.T),
                                                      reference_pressure)
"""
    DensityField(model, reference_pressure)
Return an `KernelFunctionOperation` at `(Center, Center, Center)` that computes the
potential density from the salinity and temperature tracers in `model` at `reference_pressure`.
"""
DensityField(model, reference_pressure) =
    KernelFunctionOperation{Center, Center, Center}(σ, model.grid, model,
                                                    reference_pressure)
"""
    σ_anomalyᶜᶜᶠ(i, j, k, grid, model, reference_pressure)
Compute potential density anomaly at `reference_pressure` from salinity and temperature tracers
in `model` that have been interpolated from `(Center, Center, Center)` to
`(Center, Center, Face)`.
"""
σ_anomalyᶜᶜᶠ(i, j, k, grid, model, reference_pressure, σ_reference) =
    gsw_rho(ℑzᵃᵃᶠ(i, j, k, grid, model.tracers.S), ℑzᵃᵃᶠ(i, j, k, grid, model.tracers.T),
            reference_pressure) - σ_reference
"""
    function InterpolatedDensityAnomaly(model, reference_pressure)
Return a `KernelFunctionOperation` at `(Center, Center, Face)` to compute the potential
density anomaly from the salinity and temperature tracers from `model` that have been
interpolated from `(Center, Center, Center)` to `(Center, Center, Face)`
"""
function InterpolatedDensityAnomaly(model, reference_pressure)

    σ_reference = model.buoyancy.model.equation_of_state.reference_density
    return KernelFunctionOperation{Center, Center, Face}(σ_anomalyᶜᶜᶠ,
                                                         model.grid, model,
                                                         reference_pressure,
                                                         σ_reference)

end
