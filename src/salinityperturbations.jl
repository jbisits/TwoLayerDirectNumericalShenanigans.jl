"`show` for `SalinityPerturbation`"
function Base.show(io::IO, sp::SalinityPerturbation)
    if sp isa GaussianProfile
        println(io, "$(typeof(sp))")
        println(io, " ┣━━ interface_location: z = $(sp.interface_location)m")
        println(io, " ┣━━ centre_of_Gaussian: z = $(sp.μ)m")
        println(io, " ┣━━━ width_of_Gaussian: $(sp.σ)")
        print(io,   " ┗━━━━━━━━━━━━━━━ scale: $(sp.scale)")
    elseif sp isa GaussianBlob
        println(io, "$(typeof(sp))")
        println(io, " ┣━━━━━━ depth_location: z = $(round(sp.depth; digits = 3))m")
        println(io, " ┣━━ centre_of_Gaussian: $(sp.μ)")
        println(io, " ┣━━━ width_of_Gaussian: $(sp.σ)")
        print(io,   " ┗━━━━━━━━━━━━━━━ scale: $(sp.scale)")
    elseif sp isa RandomPerturbations
        println(io, "$(typeof(sp))")
        println(io, " ┣━━ noise_location: z = $(round(sp.depth; digits = 3))m")
        print(io,   " ┗━━━━━ noise_scale: $(sp.scale)")
    end
end
"""
    struct GaussianProfile
Container for a Gaussian profile salinity perturbation in the upper layer.
"""
struct GaussianProfile{T} <: SalinityPerturbation
    "Location of interface betweenn upper and lower layers."
    interface_location :: T
    "Center of Gaussin profile in the upper layer."
    μ :: T
    "With of Gaussian profile in the upper layer."
    σ :: T
    "Scale the Gausssian profile, defaults to 1 i.e. it is a pdf."
    scale :: T
end
GaussianProfile(interface_location, μ, σ; scale = 1.0) =
    GaussianProfile(interface_location, μ, σ, scale)
"""
    struct GaussianBlob
Container for a horizontal Gaussian blob salinity perturbation at `depth` in the upper layer.
"""
struct GaussianBlob{T} <: SalinityPerturbation
    "Depth at which to set the horizontal Gaussian blob of salinity."
    depth :: T
    "Centre of the blob."
    μ :: Vector{T}
    "Width of the blob."
    σ :: T
    "Scale the Gausssian blob, defaults to 1 i.e. it is a pdf."
    scale :: T
end
GaussianBlob(interface_location, μ, σ; scale = 1.0) =
    GaussianBlob(interface_location, μ, σ, scale)
"""
    function perturb_salinity(z, salinity_perturbation::GaussianProfile)
Perturb salinity by setting a vertical Gaussian profile in the upper layer centred at `μ`
with width `σ`.
"""
function perturb_salinity(z, salinity_perturbation::GaussianProfile)

    μ, σ, scale = salinity_perturbation.μ, salinity_perturbation.σ,
                  salinity_perturbation.scale

    if z > salinity_perturbation.interface_location
        scale * exp(-(z - μ)^2 / 2*(σ)^2) / sqrt(2*π*σ^2)
    else
        0
    end

end
"""
    function perturb_salinity(x, y, z, salinity_perturbation::GaussianBlob)
Perturb salinity by setting a horizontal Gaussian blob in the upper layer centred at `μ`
with width `σ` at depth `depth`. **Note** the `depth` needs to be an exact match to a
depth at that the `Center` in the `z` direction.
"""
function perturb_salinity(x, y, z, salinity_perturbation::GaussianBlob)

    μ, σ, scale = salinity_perturbation.μ, salinity_perturbation.σ,
                  salinity_perturbation.scale

    if z == salinity_perturbation.depth
        scale * exp(- ((x - μ[1])^2 + (y - μ[2])^2) / 2*σ^2) / (2*π*σ^2)
    else
        0
    end

end
