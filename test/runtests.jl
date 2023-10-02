using TwoLayerDirectNumericalShenanigans, Test
using TwoLayerDirectNumericalShenanigans: perturb_tracer

include("twolayercontainers_test.jl")

@testset "TwoLayerDNS upper layer containers" begin
    for (i, uc) ∈ enumerate(upper_layer_containers)
        if uc isa IsohalineUpperLayerInitialConditions
            @test isequal(uc.S, S₀ᵘ[i])
            @test isequal(uc.T₀ᵘ, T₀_array[i])
        elseif uc isa IsothermalUpperLayerInitialConditions
            @test isequal(uc.S₀ᵘ, S₀ᵘ[i])
            @test isequal(uc.T, T₀_array[i])
        else
            @test isequal(uc.S₀ᵘ, S₀ᵘ[i])
            @test isequal(uc.T₀ᵘ, T₀_array[i])
        end
    end
end
@testset "TwoLayerDNS two layer containers default values" begin
    for (i, tc) ∈ enumerate(two_layer_containers)
        @test isequal(tc.S₀ᵘ, S₀ᵘ[i])
        @test isequal(tc.S₀ˡ, S₀ˡ)
        @test isequal(tc.ΔS₀, ΔS[i])
        @test isequal(tc.T₀ᵘ, T₀_array[i])
        @test isequal(tc.T₀ˡ, T₀ˡ)
        @test isequal(tc.ΔT₀, ΔT[i])
    end
end

@testset "TwoLayerDNS two layer containers non default values" begin
    for (i, tc) ∈ enumerate(two_layer_containers_alt)
        @test isequal(tc.S₀ᵘ, S₀ᵘ[i])
        @test isequal(tc.S₀ˡ, S₀ˡ_alt)
        @test isequal(tc.ΔS₀, ΔS_alt[i])
        @test isequal(tc.T₀ᵘ, T₀_array[i])
        @test isequal(tc.T₀ˡ, T₀ˡ_alt)
        @test isequal(tc.ΔT₀, ΔT_alt[i])
    end
    @test isequal(isohaline_twolayer_alt.S₀ᵘ, S₀ˡ_alt)
    @test isequal(isohaline_twolayer_alt.S₀ˡ, S₀ˡ_alt)
    @test isequal(isohaline_twolayer_alt.ΔS₀, 0)
    @test isequal(isohaline_twolayer_alt.T₀ᵘ, T₀_array[4])
    @test isequal(isohaline_twolayer_alt.T₀ˡ, T₀ˡ_alt)
    @test isequal(isohaline_twolayer_alt.ΔT₀, ΔT_alt[4])
    @test isequal(isothermal_twolayer_alt.S₀ᵘ, S₀ᵘ[5])
    @test isequal(isothermal_twolayer_alt.S₀ˡ, S₀ˡ_alt)
    @test isequal(isothermal_twolayer_alt.ΔS₀, ΔS_alt[5])
    @test isequal(isothermal_twolayer_alt.T₀ᵘ, T₀ˡ_alt)
    @test isequal(isothermal_twolayer_alt.T₀ˡ, T₀ˡ_alt)
    @test isequal(isothermal_twolayer_alt.ΔT₀, 0)
    @test TwoLayerInitialConditions(S₀ᵘ.stable, T₀ᵘ, S₀ˡ, T₀ˡ) isa StableTwoLayerInitialConditions
    @test TwoLayerInitialConditions(S₀ᵘ.unstable, T₀ᵘ, S₀ˡ, T₀ˡ) isa UnstableTwoLayerInitialConditions
end

include("initialconditions_test.jl")
@testset "Tracer Gaussian profile" begin

    for tb ∈ tracer_profile_perturbations

        model = DNS(architecture, DOMAIN_EXTENT, resolution, diffusivities;
                    reference_density = REFERENCE_DENSITY)
        dns = TwoLayerDNS(model, profile_function, initial_conditions, tracer_perturbation = tb)
        set_two_layer_initial_conditions!(dns)

        if dns.tracer_perturbation isa SalinityGaussianProfile
            @test isequal((true, trues(length(z))), tracer_profile(dns))
        elseif dns.tracer_perturbation isa TemperatureGaussianProfile
            @test isequal((trues(length(z)), true), tracer_profile(dns))
        end
    end

end
@testset "Tracer Gaussian blob" begin

    for tb ∈ tracer_blob_perturbations

        model = DNS(architecture, DOMAIN_EXTENT, resolution, diffusivities;
                    reference_density = REFERENCE_DENSITY)
        dns = TwoLayerDNS(model, profile_function, initial_conditions, tracer_perturbation = tb)
        set_two_layer_initial_conditions!(dns)

        @test isequal((true, true), tracer_blob(dns))

    end

end
@testset "Tracer noise" begin

    for tn ∈ tracer_noise_perturbations

        model = DNS(architecture, DOMAIN_EXTENT, resolution, diffusivities;
                    reference_density = REFERENCE_DENSITY)
        dns = TwoLayerDNS(model, profile_function, initial_conditions, initial_noise = tn)
        set_two_layer_initial_conditions!(dns)

        @test isequal((true, true), tracer_noise(dns))

    end
    for tnv ∈ tracer_noise_perturbations_vec

        model = DNS(architecture, DOMAIN_EXTENT, resolution, diffusivities;
                    reference_density = REFERENCE_DENSITY)
        dns = TwoLayerDNS(model, profile_function, initial_conditions, initial_noise = tnv)
        set_two_layer_initial_conditions!(dns)

        @test isequal((true, true), tracer_noise(dns))

    end

end

@testset "Step change" begin

    model = DNS(architecture, DOMAIN_EXTENT, resolution, diffusivities;
                reference_density = REFERENCE_DENSITY)
    profile_function = StepChange(z[depth_idx])
    initial_conditions = TwoLayerInitialConditions(34.551, -1.5, 34.7, 0.5)
    dns = TwoLayerDNS(model, profile_function, initial_conditions)
    set_two_layer_initial_conditions!(dns)

    @test isequal((true, true, true, true), tracer_stepchange(dns))
end

@testset "Find depth" begin
    test_depth = rand(z)
    @test isequal(test_depth, find_depth(model, test_depth))
    test_depth_range = vcat(z[10], z[20])
    @test isequal(z[10:20], find_depth(model, test_depth_range))
end

include("kernelfunction_test.jl")

@testset "Density field computation" begin
    d_idx = findfirst(z .== transition_depth)
    σᵘ = gsw_rho.(interior(model.tracers.S, rand(1:10), rand(1:10), 1:d_idx),
                  interior(model.tracers.T, rand(1:10), rand(1:10), 1:d_idx),
                  reference_pressure)
    σˡ = gsw_rho.(interior(model.tracers.S, rand(1:10), rand(1:10), d_idx+1:100),
                  interior(model.tracers.T, rand(1:10), rand(1:10), d_idx+1:100),
                  reference_pressure)
    σ_computed_profile = vcat(σᵘ, σˡ)
    @test isequal(trues(length(z)), test_density_profile(σ_field, σ_computed_profile))
end

@testset "Interpolated density anomaly" begin
    d_idx = findfirst(z .== transition_depth)
    S_interpolation = Field(KernelFunctionOperation{Center, Center, Face}(ℑzᵃᵃᶠ, model.grid,
                                                                          model.tracers.S))
    compute!(S_interpolation)
    T_interpolation = Field(KernelFunctionOperation{Center, Center, Face}(ℑzᵃᵃᶠ, model.grid,
                                                                          model.tracers.T))
    compute!(T_interpolation)
    σ_reference = dns.model.buoyancy.model.equation_of_state.reference_density
    σᵘ = gsw_rho.(interior(S_interpolation, rand(1:10), rand(1:10), 1:d_idx),
                  interior(T_interpolation, rand(1:10), rand(1:10), 1:d_idx),
                  reference_pressure) .- σ_reference
    σˡ = gsw_rho.(interior(S_interpolation, rand(1:10), rand(1:10), d_idx+1:101),
                  interior(T_interpolation, rand(1:10), rand(1:10), d_idx+1:101),
                  reference_pressure) .- σ_reference
    σ_anom_computed_profile = vcat(σᵘ, σˡ)
    @test isequal(trues(length(σ_anom_computed_profile)),
                  test_density_anom_profile(σ_anomaly_interp_field, σ_anom_computed_profile))
end
