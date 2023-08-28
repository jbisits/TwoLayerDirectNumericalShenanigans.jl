using DirectNumericalCabbelingShenanigans
using DirectNumericalCabbelingShenanigans.TwoLayerDNS
using Test

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
