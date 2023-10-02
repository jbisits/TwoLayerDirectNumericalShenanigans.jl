T₀ᵘ = -1.5
S₀ᵘ = (stable = 34.551, cabbeling = 34.568, unstable = 34.59, isohaline = 34.7, isothermal = 34.5)
T₀ˡ = 0.5
S₀ˡ = 34.7
T₀ˡ_alt = 1.0
S₀ˡ_alt = 34.8

# Stable
stable_upper = StableUpperLayerInitialConditions(S₀ᵘ.stable, T₀ᵘ)
stable_twolayer = TwoLayerInitialConditions(stable_upper)
stable_twolayer_alt = TwoLayerInitialConditions(stable_upper; T₀ˡ = T₀ˡ_alt, S₀ˡ = S₀ˡ_alt)

# Cabbeling
cabbeling_upper = CabbelingUpperLayerInitialConditions(S₀ᵘ.cabbeling, T₀ᵘ)
cabbeling_twolayer = TwoLayerInitialConditions(cabbeling_upper)
cabbeling_twolayer_alt = TwoLayerInitialConditions(cabbeling_upper; T₀ˡ = T₀ˡ_alt, S₀ˡ = S₀ˡ_alt)

# Unstable
unstable_upper = UnstableUpperLayerInitialConditions(S₀ᵘ.unstable, T₀ᵘ)
unstable_twolayer = TwoLayerInitialConditions(unstable_upper)
unstable_twolayer_alt = TwoLayerInitialConditions(unstable_upper; T₀ˡ = T₀ˡ_alt, S₀ˡ = S₀ˡ_alt)

# Isohaline
isohaline_upper = IsohalineUpperLayerInitialConditions(S₀ᵘ.isohaline, T₀ᵘ)
isohaline_twolayer = TwoLayerInitialConditions(isohaline_upper)
isohaline_upper_alt = IsohalineUpperLayerInitialConditions(S₀ˡ_alt, T₀ᵘ)
isohaline_twolayer_alt = TwoLayerInitialConditions(isohaline_upper_alt; T₀ˡ = T₀ˡ_alt)

# Isothermal
isothermal_upper = IsothermalUpperLayerInitialConditions(S₀ᵘ.isothermal, 0.5)
isothermal_twolayer = TwoLayerInitialConditions(isothermal_upper)
isothermal_upper_alt = IsothermalUpperLayerInitialConditions(S₀ᵘ.isothermal, T₀ˡ_alt)
isothermal_twolayer_alt = TwoLayerInitialConditions(isothermal_upper_alt; S₀ˡ = S₀ˡ_alt)

upper_layer_containers = (stable_upper, cabbeling_upper, unstable_upper, isohaline_upper,
                          isothermal_upper)
two_layer_containers = (stable_twolayer, cabbeling_twolayer, unstable_twolayer,
                        isohaline_twolayer, isothermal_twolayer)
two_layer_containers_alt = (stable_twolayer_alt, cabbeling_twolayer_alt, unstable_twolayer_alt)

# True values
ΔS = Array{Float64}(undef, length(S₀ᵘ))
ΔS_alt = similar(ΔS)
for (i, S) ∈ enumerate(S₀ᵘ)
    ΔS[i] = S - S₀ˡ
    ΔS_alt[i] = S - S₀ˡ_alt
end
T₀_array = vcat(fill(T₀ᵘ, 4), 0.5)
ΔT = Array{Float64}(undef, length(S₀ᵘ))
ΔT_alt = similar(ΔT)
for (i, T) ∈ enumerate(T₀_array)
    ΔT[i] = T - T₀ˡ
    ΔT_alt[i] = T - T₀ˡ_alt
end
