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
abstract type TracerNoise{D, T} <: AbstractNoise end
"""
    struct SalinityNoise
Container for adding `scale`d random noise to the salinity field at `depth`.
"""
struct SalinityNoise{D, T} <: TracerNoise{D, T}
    "Depth at which to set random salinity perturbations."
    depth :: D
    "Scale for the random noise."
    scale :: T
end
"""
    struct TemperatuerNoise
Container for adding `scale`d random noise to the salinity field at `depth`.
"""
struct TemperatureNoise{D, T} <: TracerNoise{D, T}
    "Depth at which to set random salinity perturbations."
    depth :: D
    "Scale for the random noise."
    scale :: T
end
"""
    function perturb_tracer(z, tracer_perturbation::TracerNoise)
Perturb tracer by adding `scale`d random noise to the to the `tracer` field at `depth`.
**Note** the `depth` needs to be an exact match to a depth at that the `Center` in the `z`
direction.
Either a single `depth` and `scale` or a `Vector` of `depths` and `scales` can be passed in.
"""
function perturb_tracer(z, tracer_perturbation::TracerNoise{<:Number, <:Number})

    depth, scale = tracer_perturbation
    if z == depth
        scale * randn()
    else
        0
    end

end
function perturb_tracer(z, tracer_perturbation::TracerNoise{<:AbstractVector, <:AbstractVector})

    depths, scales = tracer_perturbation
    match_depth = z .== depths
    if sum(match_depth) != 0
        reduce(+, scales[match_depth]) * randn()
    else
        0
    end

end
function perturb_tracer(z, tracer_perturbation::TracerNoise{<:CuArray, <:CuArray})

    depths, scales = tracer_perturbation
    depths = Array(depths) # move to CPU
    scales = Array(scales) # move to CPU
    match_depth = z .== depths
    if sum(match_depth) != 0
        reduce(+, scales[match_depth]) * randn()
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
