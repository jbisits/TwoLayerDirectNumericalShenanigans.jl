"`show` for `AbstractNoise`"
function Base.show(io::IO, noise::AbstractNoise)
    println(io, "$(typeof(noise))")
    println(io, "┣━━━ noise_depth: z = $(noise.depth)m ")
    print(io,   "┗━━━ noise_scale: $(noise.scale)")
end
"`iterate` for `AbstractNoise`"
Base.iterate(noise::AbstractNoise, state = 1) =
    state > length(fieldnames(typeof(noise))) ? nothing :
                                            (getfield(noise, state), state + 1)
"Tracer Noise type"
abstract type TracerNoise <: AbstractNoise end
"""
    struct SalinityNoise
Container for adding `scale`d random noise to the salinity field at `depth`.
"""
struct SalinityNoise{T} <: TracerNoise
    "Depth at which to set random salinity perturbations."
    depth :: T
    "Scale for the random noise."
    scale :: T
end
"""
    struct TemperatuerNoise
Container for adding `scale`d random noise to the salinity field at `depth`.
"""
struct TemperatureNoise{T} <: TracerNoise
    "Depth at which to set random salinity perturbations."
    depth :: T
    "Scale for the random noise."
    scale :: T
end
"""
    function perturb_salinity(z, tracer_perturbation::TracerNoise)
Perturb tracer by adding `scale`d random noise to the salinity at `depth`.
**Note** the `depth` needs to be an exact match to a depth at that the `Center` in the `z`
direction.
"""
function perturb_tracer(z, tracer_perturbation::TracerNoise)

    if z == tracer_perturbation.depth
        tracer_perturbation.scale * randn()
    else
        0
    end

end
"""
    struct VelocityNoise
Container for adding `scale`d random noise to the velocity fields at `depth`.
"""
struct VelocityNoise{T} <: AbstractNoise
    "Depth at which to set random velocity perturbations."
    depth :: T
    "Scale for the random noise."
    scale :: T
end
"""
    function perturb_velocity!(dns::TwoLayerDNS, velocity_noise::VelocityNoise;
                               horizontal = false)
Add standard normally distributed random noise given by `velocity_noise`.To only add
horizontal random noise (i.e. in the `u` and `v` velocity fields) set `true` for the
`horizontal` keyword argument.
"""
function perturb_velocity!(model::Oceananigans.AbstractModel, velocity_noise::VelocityNoise;
                           horizontal = false)

    depth, scale = velocity_noise.depth, velocity_noise.scale
    add_noise(x, y, z) = z == depth ? scale * randn() : 0

    horizontal == true ? set!(model, u = add_noise, v = add_noise) :
                         set!(model, u = add_noise, v = add_noise, w = add_noise)

    return nothing

end
