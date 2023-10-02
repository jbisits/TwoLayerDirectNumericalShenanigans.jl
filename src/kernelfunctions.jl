# Can I use a view here?
# σ_kernel_function(i, j, k, grid, S, T, density_reference_pressure) =
#     gsw_rho(S.data[i, j, k], T.data[i, j, k], density_reference_pressure)
# σ_anomaly_kernel_function(i, j, k, grid, Sᶜᶜᶠ, Tᶜᶜᶠ, density_reference_pressure, model) =
#     gsw_rho(Sᶜᶜᶠ.data[i, j, k], Tᶜᶜᶠ.data[i, j, k], density_reference_pressure) -
#     model.buoyancy.model.equation_of_state.reference_density
# σ_op = KernelFunctionOperation{Center, Center, Center}(σ_kernel_function,
#                                                    model.grid, S, T,
#                                                    density_reference_pressure)
# # Combine this into one KFO for a tracer rather than seperatre S and T.
# # Then try and pass the interpolated S and T data straight to `σ_anomaly_kernel_function`
# # rather than computing fields. It would be better if there was one function for the density
# # computation rather than needing two.
# S_interpolation = KernelFunctionOperation{Center, Center, Face}(Oceananigans.Operators.ℑzᵃᵃᶠ, model.grid, S)
# T_interpolation = KernelFunctionOperation{Center, Center, Face}(Oceananigans.Operators.ℑzᵃᵃᶠ, model.grid, T)
# σ_opᶜᶜᶠ = KernelFunctionOperation{Center, Center, Face}(σ_anomaly_kernel_function,
#                                                         model.grid, Sᶜᶜᶠ, Tᶜᶜᶠ,
#                                                         density_reference_pressure,
#                                                         model)
# ρ(i, j, k, reference_pressure) = gsw_rho(S[i, j, k], T[i, j, k], reference_pressure)
# function interpolated_tracers(i, j, k, grid, C)

#     return

# end
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
