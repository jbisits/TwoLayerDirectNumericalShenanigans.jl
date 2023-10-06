using TwoLayerDirectNumericalShenanigans: PotentialDensityField, DensityField
using Oceananigans: Operators.ℑzᵃᵃᶠ
using Oceananigans: BuoyancyModels.ρ′, BuoyancyModels.θ_and_sᴬ
using GibbsSeaWater

diffusivities = (ν = 1e-4, κ = (S = 1e-5, T = 1e-5))
resolution = (Nx = 10, Ny = 10, Nz = 100)

model = DNS(architecture, DOMAIN_EXTENT, resolution, diffusivities;
            reference_density = REFERENCE_DENSITY)

z = znodes(model.grid, Center())

T₀ᵘ = -1.5
S₀ᵘ = 34.568
cabbeling = StableUpperLayerInitialConditions(S₀ᵘ, T₀ᵘ)
initial_conditions = TwoLayerInitialConditions(cabbeling)
transition_depth = find_depth(model, INTERFACE_LOCATION)
profile_function = StepChange(transition_depth)

dns = TwoLayerDNS(model, profile_function, initial_conditions)

set_two_layer_initial_conditions!(dns)

reference_pressure = 0
parameters = (pᵣ = reference_pressure,)
model = dns.model
pd_field = PotentialDensityField(model, parameters)
compute!(pd_field)

function test_potential_density_profile(pd_field, computed_density_profile, atol)

    vertical_slice = interior(pd_field, rand(1:10), rand(1:10), :)

    return isapprox.(vertical_slice, computed_density_profile; atol)

end

d_field = DensityField(model)
compute!(d_field)

function test_density_profile(d_field, computed_density_profile, atol)

    vertical_slice = interior(d_field, rand(1:10), rand(1:10), :)

    return isapprox.(vertical_slice, computed_density_profile; atol)

end
