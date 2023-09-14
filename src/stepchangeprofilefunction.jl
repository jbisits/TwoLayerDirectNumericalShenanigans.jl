abstract type AbstractStepChangeProfileFunction <: AbstractProfileFunction end
"`show` for `AbstractStepChangeProfileFunction`"
function Base.show(io::IO, stpf::AbstractStepChangeProfileFunction)
    println(io, "$(typeof(stpf))")
    print(io, "┗━━ interface_location: z = $(stpf.interface_location)m")
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
