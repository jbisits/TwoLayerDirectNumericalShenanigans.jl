T₀ᵘ = -1.5
S₀ᵘ = (stable = 34.551, cabbeling = 34.568, unstable = 34.59, isohaline = 34.7, isothermal = 34.5)
T₀ˡ = 0.5
S₀ˡ = 34.7

# Stable
stable_upper = StableUpperLayerInitialConditions(S₀ᵘ.stable, T₀ᵘ)
stable_twolayer = TwoLayerInitialConditions(stable_upper)

# Cabbeling
cabbeling_upper = CabbelingUpperLayerInitialConditions(S₀ᵘ.cabbeling, T₀ᵘ)
cabbeling_twolayer = TwoLayerInitialConditions(cabbeling_upper)

# Unstable
unstable_upper = UnstableUpperLayerInitialConditions(S₀ᵘ.unstable, T₀ᵘ)
unstable_twolayer = TwoLayerInitialConditions(unstable_upper)

# Isohaline
isohaline_upper = IsohalineUpperLayerInitialConditions(S₀ᵘ.isohaline, T₀ᵘ)
isohaline_twolayer = TwoLayerInitialConditions(isohaline_upper)

# Isothermal
isothermal_upper = IsothermalUpperLayerInitialConditions(S₀ᵘ.isothermal, 0.5)
isothermal_twolayer = TwoLayerInitialConditions(isothermal_upper)

upper_layer_containers = (stable_upper, cabbeling_upper, unstable_upper, isohaline_upper,
                          isothermal_upper)
two_layer_containers = (stable_twolayer, cabbeling_twolayer, unstable_twolayer,
                        isohaline_twolayer, isothermal_twolayer)

# True values
ΔS = Array{Float64}(undef, length(S₀ᵘ))
for (i, S) ∈ enumerate(S₀ᵘ)
    ΔS[i] = S - S₀ˡ
end
T₀_array = vcat(fill(T₀ᵘ, 4), 0.5)
ΔT = Array{Float64}(undef, length(S₀ᵘ))
for (i, T) ∈ enumerate(T₀_array)
    ΔT[i] = T - T₀ˡ
end
