using TwoLayerDirectNumericalShenanigans, Test
using TwoLayerDirectNumericalShenanigans: perturb_tracer

include("twolayercontainers.jl")

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
end

include("initialconditions.jl")
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

end
