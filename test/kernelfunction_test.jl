diffusivities = (ν = 1e-4, κ = (S = 1e-5, T = 1e-5))
resolution = (Nx = 10, Ny = 10, Nz = 100)

model = DNSModel(architecture, DOMAIN_EXTENT, resolution, diffusivities;
                reference_density = REFERENCE_DENSITY)

T₀ᵘ = -1.5
S₀ᵘ = 34.568
cabbeling = StableUpperLayerInitialConditions(S₀ᵘ, T₀ᵘ)
initial_conditions = TwoLayerInitialConditions(cabbeling)
transition_depth = find_depth(model, INTERFACE_LOCATION)
profile_function = StepChange(transition_depth)

tldns = TwoLayerDNS(model, profile_function, initial_conditions)

set_two_layer_initial_conditions!(tldns)
