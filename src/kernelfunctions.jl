"`(Center, Center, Center)` vertical velocity `Field`"
wᶜᶜᶜ(model) = wᶜᶜᶜ(model.veolcities.w, model.grid)
wᶜᶜᶜ(w, grid) = KernelFunctionOperation{Center, Center, Center}(ℑzᵃᵃᶜ, grid, w)
"`(Center, Center, Face)` vertical buoyancy gradient `Field`"
∂b∂z(model) = ∂b∂z(model.buoyancy, model.grid, model.tracers)
∂b∂z(b, grid, tracers) = KernelFunctionOperation{Center, Center, Face}(∂z_b, grid, b, tracers)

@inline vertical_buoyancy_flux(i, j, k, grid, b::SeawaterBuoyancy, C, w) =
        -ℑzᵃᵃᶜ(i, j, k, w) * buoyancy_perturbationᶜᶜᶜ(i, j, k, grid, b, C)
@inline function vertical_buoyancy_flux(model)

    grid = model.grid
    b = model.buoyancy.model
    C = model.tracers
    w = model.velocities.w

    return KernelFunctionOperation{Center, Center, Center}(vertical_buoyancy_flux, grid, b, C, w)
end
@inline function Kᵥ(i, j, k, grid, b::SeawaterBuoyancy, C, w)

    if ∂z_b(i, j, k, grid, b, C) != 0
        (-ℑzᵃᵃᶜ(i, j, k, w) * buoyancy_perturbationᶜᶜᶜ(i, j, k, grid, b, C)) / ∂z_b(i, j, k, grid, b, C)
    else
        0
    end

end
function InferredVerticalDiffusivity(model)

    grid = model.grid
    b = model.buoyancy.model
    C = model.tracers
    w = model.velocities.w

    return KernelFunctionOperation{Center, Center, Center}(Kᵥ, grid, b, C, w)
end

## This has been implemented (by me) in Oceananigans.jl as of v0.89.3. Once I know
# everything is working I will remove this in favour of Oceananigans.jl version.

# "Extend `ρ′` to compute at user defined reference geopotential height"
# SeawaterPolynomials.ρ(i, j, k, grid, eos, θ, sᴬ) = ρ(θ_and_sᴬ(i, j, k, θ, sᴬ)..., Zᶜᶜᶜ(i, j, k, grid),eos)
# SeawaterPolynomials.ρ(i, j, k, grid, eos, θ, sᴬ, Zᵣ) = ρ(θ_and_sᴬ(i, j, k, θ, sᴬ)..., Zᵣ, eos)
# """
#     function densityᶜᶜᶜ(i, j, k, grid, b::SeawaterBuoyancy, C)
# Compute the density of seawater at grid point `(i, j, k)` using `SeawaterBuoyancy`.
# """
# @inline function densityᶜᶜᶜ(i, j, k, grid, b::SeawaterBuoyancy, C)
#     T, S = get_temperature_and_salinity(b, C)
#     return ρ(i, j, k, grid, b.equation_of_state, T, S)
# end
# density(model) = density(model.buoyancy, model.grid, model.tracers)
# density(b, grid, tracers) =
#     KernelFunctionOperation{Center, Center, Center}(densityᶜᶜᶜ, grid, b.model, tracers)
# Density(model) = density(model)
# """
#     function potential_densityᶜᶜᶜ(i, j, k, grid, b::SeawaterBuoyancy, C, parameters)
# Compute the potential density of seawater at grid point `(i, j, k)`
# at reference pressure `parameters.pᵣ` using `SeawaterBuoyancy`.
# """
# @inline function potential_densityᶜᶜᶜ(i, j, k, grid, b::SeawaterBuoyancy, C, parameters)
#     T, S = get_temperature_and_salinity(b, C)
#     Zᵣ = parameters.Zᵣ
#     return ρ(i, j, k, grid, b.equation_of_state, T, S, Zᵣ)
# end
# potential_density(model, parameters) = potential_density(model.buoyancy, model.grid,
#                                                          model.tracers, parameters)
# potential_density(b, grid, tracers, parameters) =
#     KernelFunctionOperation{Center, Center, Center}(potential_densityᶜᶜᶜ, grid, b.model,
#                                                     tracers, parameters)
# PotentialDensity(model, parameters) = potential_density(model, parameters)
