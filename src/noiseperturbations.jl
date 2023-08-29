"""
    struct SalinityNoise
Container for adding `scale`d random noise to the salinity field at `depth`.
"""
struct SalinityNoise{T} <: AbstractNoise
    "Depth at which to set random salinity perturbations."
    depth :: T
    "Scale for the random noise."
    scale :: T
end
"""
    function perturb_salinity(z, salinity_perturbation::SalinityNoise)
Perturb salinity by adding `scale`d random noise to the salinity at `depth`.
**Note** the `depth` needs to be an exact match to a depth at that the `Center` in the `z`
direction.
"""
function perturb_salinity(z, salinity_perturbation::SalinityNoise)

    if z == salinity_perturbation.depth
        salinity_perturbation.scale * randn()
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
