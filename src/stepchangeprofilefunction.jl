abstract type AbstractStepChangeProfileFunction <: AbstractProfileFunction end
"`show` for `AbstractStepChangeProfileFunction`"
function Base.show(io::IO, stpf::AbstractStepChangeProfileFunction)
    if stpf isa StepChange
        println(io, "$(typeof(stpf))")
        print(io, "┗━━ interface_location: z = $(stpf.interface_location)m")
    elseif stpf isa StepChangeLinearGradient
        println(io, "$(typeof(stpf))")
        println(io, "┣━ interface_location: z = $(stpf.interface_location)m ")
        println(io, "┣━━━━━━━━━━━━━━━ dSdz: $(stpf.dSdz) ")
          print(io, "┗━━━━━━━━━━━━━━━ dTdZ: $(stpf.dTdZ)")
    end
end
"`iterate` for `AbstractStepChangeProfileFunction`"
Base.iterate(pf::AbstractStepChangeProfileFunction, state = 1) =
    state > length(fieldnames(typeof(pf))) ? nothing :
                                            (getfield(pf, state), state + 1)
"""
    struct StepChange
Containter with the depth of the `interface_location` between the upper and lower layers.
The resulting profile is a step change between the upper and lower layer salinity and
temperature values.
"""
struct StepChange{T} <: AbstractStepChangeProfileFunction
    "Location of the interface between the two layers"
    interface_location :: T
end
Heaviside(z, Cˡ::Number, Cᵘ::Number, profile_function::StepChange) =
    z > profile_function.interface_location ? Cᵘ : Cˡ

"""
    struct StepChangeLinearGradient
Containter with the depth of the `interface_location` between the upper and lower layers
and values to construct a linear gradient for salinity (`dSdz`) and temperature (`dTdz`).
The resulting profile is a step change between the upper and lower layer salinity and
temperature that has a linear (with depth) gradient in each layer.
"""
struct StepChangeLinearGradient{T} <: AbstractStepChangeProfileFunction
    "Location of the interface between the two layers"
    interface_location :: T
    "Linear gradient value for salinity"
                  dSdz :: T
    "Salinity offset to ensure salinity values at interface are S₀ᵘ and S₀ˡ"
       salinity_offset :: Tuple{T, T}
    "Linear gradient value for temperature"
                  dTdZ :: T
    "Temperature offset to ensure salinity values at interface are T₀ᵘ and T₀ˡ"
    temperature_offset :: Tuple{T, T}
end
function StepChangeLinearGradient(interface_location, dSdz, dTdz, model)

    offset_depth = model.architecture isa CPU ?  begin
                                                z = znodes(model.grid, Center(), Center(), Center())
                                                il_idx = findfirst(z .> interface_location)
                                                z[il_idx], z[il_idx - 1]
                                            end :
                                                allowscalar() do
                                                z = znodes(model.grid, Center(), Center(), Center())
                                                il_idx = findfirst(z .> interface_location)
                                                z[il_idx], z[il_idx - 1]
                                            end

    S_offset = (dSdz * offset_depth[1], dSdz * offset_depth[2])
    T_offset = (dTdz * offset_depth[1], dTdz * offset_depth[2])

    return StepChangeLinearGradient(interface_location, dSdz, S_offset, dTdz, T_offset)

end
function Heaviside_with_linear_gradient(z, Cˡ::Number, Cᵘ::Number,
                                        profile_function::StepChangeLinearGradient;
                                        tracer = :S)
    interface_location = profile_function.interface_location
    dCdz = tracer == :S ? -profile_function.dSdz : profile_function.dTdZ
    upper_offset, lower_offset = tracer == :S ? profile_function.salinity_offset :
                                                profile_function.temperature_offset
    Cᵘ_offset = Cᵘ + upper_offset
    Cˡ_offset = Cˡ + lower_offset
    if z > interface_location
        Cᵘ_offset + dCdz * z
    else
        Cˡ_offset + dCdz * z
    end

end
