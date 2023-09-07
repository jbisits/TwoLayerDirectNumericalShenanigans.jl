"Iterate for `AbstractInitialConditions`"
Base.iterate(ics::AbstractInitialConditions, state = 1) =
    state > length(fieldnames(typeof(ics))) ? nothing :
                                            (getfield(ics, state), state + 1)
"Abstract super type for initial temperature and salinity in the upper layer."
abstract type UpperLayerInitialConditions <: AbstractInitialConditions end
"`show` for `UpperLayerInitialConditions`"
function Base.show(io::IO, ulic::UpperLayerInitialConditions)
    if ulic isa IsothermalUpperLayerInitialConditions
        println(io, "$(typeof(ulic))")
        println(io, "┣━━ S = $(ulic.S₀ᵘ)")
        print(io,   "┗━━ T = $(ulic.T)")
    elseif ulic isa IsohalineUpperLayerInitialConditions
        println(io, "$(typeof(ulic))")
        println(io, "┣━━ S = $(ulic.S)")
        print(io,   "┗━━ T = $(ulic.T₀ᵘ)")
    else
        println(io, "$(typeof(ulic))")
        println(io, "┣━━ S = $(ulic.S₀ᵘ)")
        print(io,   "┗━━ T = $(ulic.T₀ᵘ)")
    end
end
"""
    struct StableUpperLayerInitialConditions
Container for initial salinity and temperature conditions that are stable relative to `S₀ˡ`
and `T₀ˡ`.
"""
struct StableUpperLayerInitialConditions{T} <: UpperLayerInitialConditions
    "Initial salinity in the upper layer"
    S₀ᵘ :: T
    "Initial temperature in the upper layer"
    T₀ᵘ :: T
end
"""
    struct CabbelingUpperLayerInitialConditions
Container for initial salinity and temperature conditions that are unstable to cabbeling
relative to `S₀ˡ` and `T₀ˡ`.
"""
struct CabbelingUpperLayerInitialConditions{T} <: UpperLayerInitialConditions
    "Initial salinity in the upper layer"
    S₀ᵘ :: T
    "Initial temperature in the upper layer"
    T₀ᵘ :: T
end
"""
    struct UnstableUpperLayerInitialConditions
Container for initial salinity and temperature conditions that are unstable relative to `S₀ˡ`
and `T₀ˡ`.
"""
struct UnstableUpperLayerInitialConditions{T} <: UpperLayerInitialConditions
    "Initial salinity in the upper layer"
    S₀ᵘ :: T
    "Initial temperature in the upper layer"
    T₀ᵘ :: T
end
"""
    struct IsohalineUpperLayerInitialConditions
Container for isohaline initial salinity at (`S₀ˡ`) and initial temperature conditions `T₀ˡ`.
"""
struct IsohalineUpperLayerInitialConditions{T} <: UpperLayerInitialConditions
    "Initial (uniform) salinity over the domain"
    S   :: T
    "Initial temperature in the upper layer"
    T₀ᵘ :: T
end
"""
    struct IsohalineUpperLayerInitialConditions
Container for isohaline initial salinity at (`S₀ˡ`) and initial temperature conditions `T₀ˡ`.
"""
struct IsothermalUpperLayerInitialConditions{T} <: UpperLayerInitialConditions
    "Initial salinity in the upper layer"
    S₀ᵘ :: T
    "Initial (uniform) temperature over the domain"
    T   :: T
end
"Abstract super type for two layer initial conditions for the dns."
abstract type TwoLayerInitialConditions <: AbstractInitialConditions end
"`show` for `TwoLayerInitialConditions`"
function Base.show(io::IO, tlic::TwoLayerInitialConditions)
    println(io, "$(typeof(tlic))")
    println(io, "┣━━ upper_layer: S = $(tlic.S₀ᵘ), T = $(tlic.T₀ᵘ)")
    print(io,   "┗━━ lower_layer: S = $(tlic.S₀ˡ), T = $(tlic.T₀ˡ)")
end
"""
    function TwoLayerInitialConditions(S₀ᵘ, T₀ᵘ, S₀ˡ, T₀ˡ; reference_pressure = 0)
Return either `Stable` of `Unstable` `TwoLayerInitialConditions` depending on the
density differenece
```math
    Δρ = ρ \\left(S₀ᵘ, T₀ᵘ, pᵣ \\right) - ρ \\left(S₀ˡ, T₀ˡ, pᵣ \\right)
```
calculate from the salinity and temperature set for the upper and lower layer.
"""
function TwoLayerInitialConditions(S₀ᵘ, T₀ᵘ, S₀ˡ, T₀ˡ; reference_pressure = 0)

    Δρ = gsw_rho(S₀ᵘ, T₀ᵘ, reference_pressure) - gsw_rho(S₀ˡ, T₀ˡ, reference_pressure)

    tlics = Δρ ≤ 0 ? StableTwoLayerInitialConditions(S₀ᵘ, S₀ˡ, S₀ᵘ - S₀ˡ,
                                                     T₀ᵘ,  T₀ˡ, T₀ᵘ - T₀ˡ) :
                     UnstableTwoLayerInitialConditions(S₀ᵘ, S₀ˡ, S₀ᵘ - S₀ˡ,
                                                       T₀ᵘ,  T₀ˡ, T₀ᵘ - T₀ˡ)

    return tlics

end
"""
    struct StableTwoLayerInitialConditions
Container for initial salinity and temperature conditions that are stable.
"""
struct StableTwoLayerInitialConditions{T} <: TwoLayerInitialConditions
    "Initial salinity in the upper layer"
    S₀ᵘ :: T
    "Initial salinity in the lower layer"
    S₀ˡ :: T
    "Initial difference in salinity between the layers"
    ΔS₀ :: T
    "Initial temperature in the upper layer"
    T₀ᵘ :: T
    "Initial temperature in the lower layer"
    T₀ˡ :: T
    "Initial temperature difference between the layers"
    ΔT₀ :: T
end
TwoLayerInitialConditions(initial_conditions::StableUpperLayerInitialConditions;
                          S₀ˡ = 34.7, T₀ˡ = 0.5) =
    StableTwoLayerInitialConditions(initial_conditions.S₀ᵘ, S₀ˡ,
                                    initial_conditions.S₀ᵘ - S₀ˡ, initial_conditions.T₀ᵘ,
                                    T₀ˡ, initial_conditions.T₀ᵘ -T₀ˡ)
"""
    struct CabbelingTwoLayerInitialConditions
Container for initial salinity and temperature conditions that are unstable to cabbeling.
"""
struct CabbelingTwoLayerInitialConditions{T} <: TwoLayerInitialConditions
    "Initial salinity in the upper layer"
    S₀ᵘ :: T
    "Initial salinity in the lower layer"
    S₀ˡ :: T
    "Initial difference in salinity between the layers"
    ΔS₀ :: T
    "Initial temperature in the upper layer"
    T₀ᵘ :: T
    "Initial temperature in the lower layer"
    T₀ˡ :: T
    "Initial temperature difference between the layers"
    ΔT₀ :: T
end
TwoLayerInitialConditions(initial_conditions::CabbelingUpperLayerInitialConditions;
                          S₀ˡ = 34.7, T₀ˡ = 0.5) =
    CabbelingTwoLayerInitialConditions(initial_conditions.S₀ᵘ, S₀ˡ,
                                       initial_conditions.S₀ᵘ - S₀ˡ, initial_conditions.T₀ᵘ,
                                       T₀ˡ, initial_conditions.T₀ᵘ -T₀ˡ)
"""
    struct UnstableTwoLayerInitialConditions
Container for initial salinity and temperature conditions that are gravitationally unstable.
"""
struct UnstableTwoLayerInitialConditions{T} <: TwoLayerInitialConditions
    "Initial salinity in the upper layer"
    S₀ᵘ :: T
    "Initial salinity in the lower layer"
    S₀ˡ :: T
    "Initial difference in salinity between the layers"
    ΔS₀ :: T
    "Initial temperature in the upper layer"
    T₀ᵘ :: T
    "Initial temperature in the lower layer"
    T₀ˡ :: T
    "Initial temperature difference between the layers"
    ΔT₀ :: T
end
TwoLayerInitialConditions(initial_conditions::UnstableUpperLayerInitialConditions;
                          S₀ˡ = 34.7, T₀ˡ = 0.5) =
    UnstableTwoLayerInitialConditions(initial_conditions.S₀ᵘ, S₀ˡ,
                                      initial_conditions.S₀ᵘ - S₀ˡ, initial_conditions.T₀ᵘ,
                                      T₀ˡ, initial_conditions.T₀ᵘ -T₀ˡ)
"""
    struct IsohalineTwoLayerInitialConditions
Container for initial salinity and temperature conditions where salinity is uniform over
domain.
"""
struct IsohalineTwoLayerInitialConditions{T} <: TwoLayerInitialConditions
    "Initial salinity in the upper layer"
    S₀ᵘ :: T
    "Initial salinity in the lower layer"
    S₀ˡ :: T
    "Initial difference in salinity between the layers"
    ΔS₀ :: T
    "Initial temperature in the upper layer"
    T₀ᵘ :: T
    "Initial temperature in the lower layer"
    T₀ˡ :: T
    "Initial temperature difference between the layers"
    ΔT₀ :: T
end
TwoLayerInitialConditions(initial_conditions::IsohalineUpperLayerInitialConditions;
                          T₀ˡ = 0.5) =
    IsohalineTwoLayerInitialConditions(initial_conditions.S, initial_conditions.S,
                                       0.0, initial_conditions.T₀ᵘ, T₀ˡ,
                                       initial_conditions.T₀ᵘ -T₀ˡ)
                                       """
    struct IsothermalTwoLayerInitialConditions
Container for initial salinity and temperature conditions where temperature is uniform over
the domain.
"""
struct IsothermalTwoLayerInitialConditions{T} <: TwoLayerInitialConditions
    "Initial salinity in the upper layer"
    S₀ᵘ :: T
    "Initial salinity in the lower layer"
    S₀ˡ :: T
    "Initial difference in salinity between the layers"
    ΔS₀ :: T
    "Initial temperature in the upper layer"
    T₀ᵘ :: T
    "Initial temperature in the lower layer"
    T₀ˡ :: T
    "Initial temperature difference between the layers"
    ΔT₀ :: T
end
TwoLayerInitialConditions(initial_conditions::IsothermalUpperLayerInitialConditions;
                          S₀ˡ = 34.7) =
    IsothermalTwoLayerInitialConditions(initial_conditions.S₀ᵘ, S₀ˡ,
                                        initial_conditions.S₀ᵘ - S₀ˡ, initial_conditions.T,
                                        initial_conditions.T, 0.0)
