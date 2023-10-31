function run_sim(save_file)

    diffusivities = (ν = 0, κ = (S = 0, T = 0))
    resolution = (Nx = 10, Ny = 10, Nz = 100)

    model = DNSModel(architecture, DOMAIN_EXTENT, resolution, diffusivities;
                    reference_density = REFERENCE_DENSITY)

    T₀ᵘ = -1.5
    S₀ᵘ = 34.58
    cabbeling = CabbelingUpperLayerInitialConditions(S₀ᵘ, T₀ᵘ)
    initial_conditions = TwoLayerInitialConditions(cabbeling)
    transition_depth = find_depth(model, INTERFACE_LOCATION)
    td_idx = findfirst(znodes(model.grid, Center()) .== transition_depth)
    profile_function = StepChange(transition_depth)
    tldns = TwoLayerDNS(model, profile_function, initial_conditions)

    set_two_layer_initial_conditions!(tldns)

    ## build the simulation
    Δt = 1
    stop_time = 10
    save_schedule = 5 # seconds
    output_path = joinpath(@__DIR__, "outputs/")
    simulation = TLDNS_simulation_setup(tldns, Δt, stop_time, save_schedule; save_file, output_path)

    run!(simulation)

    return simulation, td_idx, tldns
end
