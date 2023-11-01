using TwoLayerDirectNumericalShenanigans, Test
using TwoLayerDirectNumericalShenanigans: perturb_tracer
using Oceananigans.Fields
using SeawaterPolynomials
using NCDatasets, JLD2
import SeawaterPolynomials.ρ
using GibbsSeaWater: gsw_p_from_z
using CUDA: has_cuda_gpu

architecture = has_cuda_gpu() ? GPU() : CPU()
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

        model = DNSModel(architecture, DOMAIN_EXTENT, resolution, diffusivities)
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

        model = DNSModel(architecture, DOMAIN_EXTENT, resolution, diffusivities)
        dns = TwoLayerDNS(model, profile_function, initial_conditions, tracer_perturbation = tb)
        set_two_layer_initial_conditions!(dns)

        @test isequal((true, true), tracer_blob(dns))

    end

end
@testset "Tracer noise" begin

    for tn ∈ tracer_noise_perturbations

        model = DNSModel(architecture, DOMAIN_EXTENT, resolution, diffusivities)
        dns = TwoLayerDNS(model, profile_function, initial_conditions, initial_noise = tn)
        set_two_layer_initial_conditions!(dns)

        @test isequal((true, true), tracer_noise(dns))

    end
    for tnv ∈ tracer_noise_perturbations_vec

        model = DNSModel(architecture, DOMAIN_EXTENT, resolution, diffusivities)
        dns = TwoLayerDNS(model, profile_function, initial_conditions, initial_noise = tnv)
        set_two_layer_initial_conditions!(dns)

        @test isequal((true, true), tracer_noise(dns))

    end

end

@testset "Step change" begin

    model = DNSModel(architecture, DOMAIN_EXTENT, resolution, diffusivities)
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

include("output_test.jl")
@testset "Saving and computed output" begin

    @testset "NetCDF" begin
        simulation, td, tldns = run_sim(:netcdf)
        NCDataset(simulation.output_writers[:tracers].filepath) do ds
            S, T = ds[:S], ds[:T]
            @test all(S[:, :, 1:td, :]     .== tldns.initial_conditions.S₀ˡ)
            @test all(S[:, :, td+1:end, :] .== tldns.initial_conditions.S₀ᵘ)
            @test all(T[:, :, 1:td, :]     .== tldns.initial_conditions.T₀ˡ)
            @test all(T[:, :, td+1:end, :] .== tldns.initial_conditions.T₀ᵘ)
        end
        NCDataset(simulation.output_writers[:computed_output].filepath) do ds
            σ = ds[:σ]
            @test all(σ[:, :, 1:td, :] .== SeawaterPolynomials.ρ(tldns.initial_conditions.T₀ˡ,
                                             tldns.initial_conditions.S₀ˡ, 0,
                                             tldns.model.buoyancy.model.equation_of_state))
            @test all(σ[:, :, td+1:end, :] .== SeawaterPolynomials.ρ(tldns.initial_conditions.T₀ᵘ,
                                                 tldns.initial_conditions.S₀ᵘ, 0,
                                                 tldns.model.buoyancy.model.equation_of_state))
        end
    end

    @testset "jld2" begin
        simulation, td, tldns = run_sim(:jld2)
        S = FieldTimeSeries(simulation.output_writers[:tracers].filepath, "S")
        T = FieldTimeSeries(simulation.output_writers[:tracers].filepath, "T")
        @test all(S.data[:, :, 1:td, :]     .== tldns.initial_conditions.S₀ˡ)
        @test all(S.data[:, :, td+1:end, :] .== tldns.initial_conditions.S₀ᵘ)
        @test all(T.data[:, :, 1:td, :]     .== tldns.initial_conditions.T₀ˡ)
        @test all(T.data[:, :, td+1:end, :] .== tldns.initial_conditions.T₀ᵘ)
        σ = FieldTimeSeries(simulation.output_writers[:computed_output].filepath, "σ")
        @test all(σ.data[:, :, 1:td, :] .==
                SeawaterPolynomials.ρ(tldns.initial_conditions.T₀ˡ,
                                      tldns.initial_conditions.S₀ˡ, 0,
                                      tldns.model.buoyancy.model.equation_of_state))
        @test all(σ.data[:, :, td+1:end, :] .==
                SeawaterPolynomials.ρ(tldns.initial_conditions.T₀ᵘ,
                                      tldns.initial_conditions.S₀ᵘ, 0,
                                      tldns.model.buoyancy.model.equation_of_state))
    end

end

include("kernelfunction_test.jl")

@testset "Vertical temperature flux" begin

    T_anom = Field(TLDNS.tracer_perturbation(tldns.model, tldns.model.tracers.T))
    compute!(T_anom)
    computed_mean = sum(interior(model.tracers.T)) / *(size(interior(model.tracers.T))...)
    computed_T_anom = interior(model.tracers.T) .- computed_mean
    @test all(computed_T_anom .== interior(T_anom))
    vtflux = Field(TLDNS.vertical_tracer_flux(tldns.model, tldns.model.tracers.T))
    compute!(vtflux)
    @test all(vtflux.data .≈ 0 )

end
