"`show` for `AbstractTracerPerturbation`"
function Base.show(io::IO, tp::AbstractTracerPerturbation)
    if tp isa SalinityGaussianProfile
        println(io, "$(typeof(tp))")
        println(io, " ┣━━ interface_location: z = $(tp.interface_location)m")
        println(io, " ┣━━ centre_of_Gaussian: z = $(tp.μ)m")
        println(io, " ┣━━━ width_of_Gaussian: $(tp.σ)")
        print(io,   " ┗━━━━━━━━━━━━━━━ scale: $(tp.scale)")
    elseif tp isa SalinityGaussianBlob
        println(io, "$(typeof(tp))")
        println(io, " ┣━━━━━━ depth_location: z = $(round(tp.depth; digits = 3))m")
        println(io, " ┣━━ centre_of_Gaussian: $(tp.μ)")
        println(io, " ┣━━━ width_of_Gaussian: $(tp.σ)")
        print(io,   " ┗━━━━━━━━━━━━━━━ scale: $(tp.scale)")
    elseif tp isa SalinityNoise
        println(io, "$(typeof(tp))")
        println(io, " ┣━━ noise_location: z = $(round(tp.depth; digits = 3))m")
        print(io,   " ┗━━━━━ noise_scale: $(tp.scale)")
    end
end
"Abstract type for salinity perturbation"
abstract type SalinityPerturbation <: AbstractTracerPerturbation end
"""
    struct SalinityGaussianProfile
Container for a Gaussian profile salinity perturbation in the upper layer.
"""
struct SalinityGaussianProfile{T} <: SalinityPerturbation
    "Location of interface betweenn upper and lower layers."
    interface_location :: T
    "Center of Gaussin profile in the upper layer."
    μ :: T
    "With of Gaussian profile in the upper layer."
    σ :: T
    "Scale the Gausssian profile, defaults to 1 i.e. it is a pdf."
    scale :: T
end
SalinityGaussianProfile(interface_location, μ, σ; scale = 1.0) =
    SalinityGaussianProfile(interface_location, μ, σ, scale)
"""
    struct SalinityGaussianBlob
Container for a horizontal Gaussian blob salinity perturbation at `depth`.
"""
struct SalinityGaussianBlob{T} <: SalinityPerturbation
    "Depth at which to set the horizontal Gaussian blob of salinity."
    depth :: T
    "Centre of the blob."
    μ :: Vector{T}
    "Width of the blob."
    σ :: T
    "Scale the Gausssian blob, defaults to 1 i.e. it is a pdf."
    scale :: T
end
SalinityGaussianBlob(interface_location, μ, σ; scale = 1.0) =
    SalinityGaussianBlob(interface_location, μ, σ, scale)

abstract type TemperaturePerturbation <: AbstractTracerPerturbation end
"""
    struct TemperatureGaussianProfile
Container for a Gaussian profile temperature perturbation in the upper layer.
"""
struct TemperatureGaussianProfile{T} <: TemperaturePerturbation
    "Location of interface betweenn upper and lower layers."
    interface_location :: T
    "Center of Gaussin profile in the upper layer."
    μ :: T
    "With of Gaussian profile in the upper layer."
    σ :: T
    "Scale the Gausssian profile, defaults to 1 i.e. it is a pdf."
    scale :: T
end
TemperatureGaussianProfile(interface_location, μ, σ; scale = 1.0) =
    TemperatureGaussianProfile(interface_location, μ, σ, scale)
"""
    struct TemperatureGaussianBlob
Container for a horizontal Gaussian blob temperature perturbation at `depth`.
"""
struct TemperatureGaussianBlob{T} <: TemperaturePerturbation
    "Depth at which to set the horizontal Gaussian blob of salinity."
    depth :: T
    "Centre of the blob."
    μ :: Vector{T}
    "Width of the blob."
    σ :: T
    "Scale the Gausssian blob, defaults to 1 i.e. it is a pdf."
    scale :: T
end
TemperatureGaussianBlob(interface_location, μ, σ; scale = 1.0) =
    TemperatureGaussianBlob(interface_location, μ, σ, scale)
"""
    function perturb_tracer(z, tracer_perturbation::SalinityGaussianProfile)
Perturb salinity by setting a vertical Gaussian profile in the upper layer centred at `μ`
with width `σ`.
"""
function perturb_tracer(z, tracer_perturbation::SalinityGaussianProfile)

    μ, σ, scale = tracer_perturbation.μ, tracer_perturbation.σ,
                  tracer_perturbation.scale

    if z > tracer_perturbation.interface_location
        scale * exp(-(z - μ)^2 / 2*(σ)^2) / sqrt(2*π*σ^2)
    else
        0
    end

end
"""
    function perturb_tracer(x, y, z, tracer_perturbation::SalinityGaussianBlob)
Perturb salinity by setting a horizontal Gaussian blob in the upper layer centred at `μ`
with width `σ` at depth `depth`. **Note** the `depth` needs to be an exact match to a
depth at that the `Center` in the `z` direction.
"""
function perturb_tracer(x, y, z, tracer_perturbation::SalinityGaussianBlob)

    μ, σ, scale = tracer_perturbation.μ, tracer_perturbation.σ,
                  tracer_perturbation.scale

    if z == tracer_perturbation.depth
        scale * exp(- ((x - μ[1])^2 + (y - μ[2])^2) / 2*σ^2) / (2*π*σ^2)
    else
        0
    end

end
