abstract type AbstractContinuousProfileFunction <: AbstractProfileFunction end
"`show` for `AbstractContinuousProfileFunction`"
function Base.show(io::IO, cpf::AbstractContinuousProfileFunction)
    if cpf isa HyperbolicTangent
        println(io, "$(typeof(cpf))")
        println(io, "┣━━━━━━━━━━ interface_location: z = $(cpf.interface_location)m ")
        print(io,   "┗━━ interface_transition_scale: $(cpf.interface_transition_width)")
    elseif cpf isa Erf
        println(io, "$(typeof(cpf))")
        println(io, "┣━━ interface_location: z = $(cpf.interface_location)m")
        print(io,   "┗━━━━━━━━ time_evolved: t = $(cpf.time)s")
    elseif cpf isa MidPoint
        println(io, "$(typeof(cpf))")
        print(io, "┗━━ interface_location: z = $(cpf.interface_location)m")
    end
end
"`iterate` for `AbstractContinuousProfileFunction`"
Base.iterate(pf::AbstractContinuousProfileFunction, state = 1) =
    state > length(fieldnames(typeof(pf))) ? nothing :
                                            (getfield(pf, state), state + 1)
"""
    struct HyperbolicTangent
Container for a hyperbolic tangent profile. The `interface_transition_width` sets the width
of the transition between the upper and lower layer.
"""
struct HyperbolicTangent{T} <: AbstractContinuousProfileFunction
    "Location of the interface between the two layers."
    interface_location :: T
    "Scale the transition between the upper and lower layer salinity and temperature."
    interface_transition_width :: T
end
"""
    struct Erf
Container for a profile that is an error function. The `time` is the time which to evaluate
`erf_tracer_solution`.
"""
struct Erf{T} <: AbstractContinuousProfileFunction
    "Location of the interface between the two layers."
    interface_location :: T
    "Time at which to evaluate the error function which is solution to 1D evolution of S or T."
    time :: T
end
"""
    struct MidPoint
Container for a profile that is a linear transition through the midpoint of the salinity
and temperature between the upper and lower layers.
"""
struct MidPoint{T} <: AbstractContinuousProfileFunction
    "Location if the interface between the two layers."
    interface_location :: T
end
"""
    function erf_tracer_solution(z, C::Number, ΔC::Number, profile_function::Erf)
Solution to the heat equation for a tracer concentration field `C` subject to initial
conditions that are a Heaviside (or modified Heaviside) step function aat time `t`.

## Function arguments:

- `z` for the Oceananigans model grid to evaulate the function at;
- `C` tracer value in deeper part of the step;
- `ΔC` difference in tracer between the steps;
- `κ` the diffusivity of the tracer;
- `profile_function::Erf` container with the `interface_location` and time `t` at which to
evaluate the error function solution

Default behaviour puts the `interface_location` in the centre of the depth range given by `z`.
"""
erf_tracer_solution(z, Cₗ::Number, ΔC::Number, κ::Number, profile_function::Erf) =
    Cₗ + 0.5 * ΔC * (1 + erf((z - profile_function.interface_location) /
                              sqrt(4 * κ * profile_function.time)))
"""
    tanh_initial_condition(z, Cᵤ::Number, ΔC::Number, interface_location)
Set a hyperbolic tangent initial condition for a tracer `C` over the vertical domain `z`.

## Function arguments

- `z` for the Oceananigans model grid to evaluate the function at;
- `Cˡ` tracer value in the lower layer;
- `ΔC` difference in tracer between upper layer and lower layer;
- `profile_function::HyperbolicTangent` container with `interface_location` and the
`interface_transition_width`.
"""
tanh_initial_condition(z, Cˡ::Number, ΔC::Number, profile_function::HyperbolicTangent) =
    Cˡ + 0.5 * ΔC * (1  + tanh(profile_function.interface_transition_width *
                               (z - profile_function.interface_location)))
"""
    midpoint(z, Cˡ::Number, Cᵘ::Number, profile_function::MidPoint)
Set the change between the upper and lower layer as the midpoint between `Cˡ` and `Cᵘ`.
This midpoint change takes place over three grid cells (in the vertical) with the midpoint
being set at `profile_location.interface_location`.

## Function arguments

- `z` for the Oceananigans model grid to evaluate the function at;
- `Cˡ` tracer value in the lower layer;
- `Cᵘ` tracer value in the upper layer;
- `profile_function::MidPoint` container with `interface_location` - that is where to set
midpoint transition.
"""
function midpoint(z, Cˡ::Number, Cᵘ::Number, profile_function::MidPoint)

    if z > profile_function.interface_location
        Cᵘ
    elseif z < profile_function.interface_location
        Cˡ
    elseif z == profile_function.interface_location
        0.5 * (Cᵘ + Cˡ)
    end

end
