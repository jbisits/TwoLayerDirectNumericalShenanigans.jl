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
@testset "TwoLayerDNS two layer containers" begin
    for (i, tc) ∈ enumerate(two_layer_containers)
        @test isequal(tc.S₀ᵘ, S₀ᵘ[i])
        @test isequal(tc.S₀ˡ, S₀ˡ)
        @test isequal(tc.ΔS₀, ΔS[i])
        @test isequal(tc.T₀ᵘ, T₀_array[i])
        @test isequal(tc.T₀ˡ, T₀ˡ)
        @test isequal(tc.ΔT₀, ΔT[i])
    end
end
