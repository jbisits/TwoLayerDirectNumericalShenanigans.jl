using TwoLayerDirectNumericalShenanigans: DensityField, InterpolatedDensityAnomaly
using Oceananigans: Operators.ℑzᵃᵃᶠ
using GibbsSeaWater

architecture = CPU()
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
σ_field = Field(DensityField(dns.model, 0))
compute!(σ_field)

function test_density_profile(σ_field, computed_density_profile)

    vertical_slice = interior(σ_field, rand(1:10), rand(1:10), :)

    return vertical_slice .== computed_density_profile

end

σ_anomaly_interp_field = Field(InterpolatedDensityAnomaly(dns.model, reference_pressure))
compute!(σ_anomaly_interp_field)

function test_density_anom_profile(σ_anomaly_interp_field, computed_density_profile)

    vertical_slice = interior(σ_anomaly_interp_field, rand(1:10), rand(1:10), :)

    return vertical_slice .== computed_density_profile

end
