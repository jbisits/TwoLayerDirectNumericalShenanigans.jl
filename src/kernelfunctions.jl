"Extend `ρ′` to compute at user defined reference pressure"
Oceananigans.BuoyancyModels.ρ′(i, j, k, grid, eos, θ, sᴬ, pᵣ) = ρ′(θ_and_sᴬ(i, j, k, θ, sᴬ)..., pᵣ, eos)
"""
    function densityᶜᶜᶜ(i, j, k, grid, b::SeawaterBuoyancy, C)
Compute the density of seawater at grid point `(i, j, k)` using `SeawaterBuoyancy`.
"""
@inline function densityᶜᶜᶜ(i, j, k, grid, b::SeawaterBuoyancy, C)
    T, S = get_temperature_and_salinity(b, C)
    return  b.equation_of_state.reference_density + ρ′(i, j, k, grid, b.equation_of_state, T, S)
end
density(model) = density(model.buoyancy, model.grid, model.tracers)
density(b, grid, tracers) =
    KernelFunctionOperation{Center, Center, Center}(densityᶜᶜᶜ, grid, b.model, tracers)
DensityField(model) = Field(density(model))
"""
    function potential_densityᶜᶜᶜ(i, j, k, grid, b::SeawaterBuoyancy, C, parameters)
Compute the potential density of seawater at grid point `(i, j, k)`
at reference pressure `parameters.pᵣ` using `SeawaterBuoyancy`.
"""
@inline function potential_densityᶜᶜᶜ(i, j, k, grid, b::SeawaterBuoyancy, C, parameters)
    T, S = get_temperature_and_salinity(b, C)
    pᵣ = parameters.pᵣ
    return  b.equation_of_state.reference_density + ρ′(i, j, k, grid, b.equation_of_state, T, S, pᵣ)
end
potential_density(model, parameters) = potential_density(model.buoyancy, model.grid,
                                                         model.tracers, parameters)
potential_density(b, grid, tracers, parameters) =
    KernelFunctionOperation{Center, Center, Center}(potential_densityᶜᶜᶜ, grid, b.model,
                                                    tracers, parameters)
PotentialDensityField(model, parameters) = Field(potential_density(model, parameters))

"`(Center, Center, Center)` vertical velocity `Field`"
wᶜᶜᶜ(model) = wᶜᶜᶜ(model.veolcities.w, model.grid)
wᶜᶜᶜ(w, grid) = KernelFunctionOperation{Center, Center, Center}(ℑzᵃᵃᶜ, grid, w)
"`(Center, Center, Face)` vertical buoyancy gradient `Field`"
∂b∂z(model) = ∂b∂z(model.buoyancy, model.grid, model.tracers)
∂b∂z(b, grid, tracers) = KernelFunctionOperation{Center, Center, Face}(∂z_b, grid, b, tracers)
∂b∂zField(model) = Field(∂b∂z(model))
