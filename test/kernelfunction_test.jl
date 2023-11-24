diffusivities = (ν = 1e-4, κ = (S = 1e-5, T = 1e-5))
resolution = (Nx = 10, Ny = 10, Nz = 100)

model = DNSModel(architecture, DOMAIN_EXTENT, resolution, diffusivities)

T₀ᵘ = -1.5
S₀ᵘ = 34.568
cabbeling = StableUpperLayerInitialConditions(S₀ᵘ, T₀ᵘ)
initial_conditions = TwoLayerInitialConditions(cabbeling)
transition_depth = find_depth(model, INTERFACE_LOCATION)
profile_function = StepChange(transition_depth)

tldns = TwoLayerDNS(model, profile_function, initial_conditions)

set_two_layer_initial_conditions!(tldns)

Eₚ_model = DNSModel(architecture, DOMAIN_EXTENT, resolution, diffusivities)
set!(Eₚ_model, T = -1.5, S = 34.568)
ρ_model = Field(seawater_density(Eₚ_model))
compute!(ρ_model)
z_grid = Field(Oceananigans.Models.model_geopotential_height(Eₚ_model))
compute!(z_grid)
dV = xspacings(model.grid, Center()) * yspacings(model.grid, Center()) * zspacings(model.grid, Center())
computed_Eₚ = model.buoyancy.model.gravitational_acceleration * ρ_model.data .* z_grid.data
