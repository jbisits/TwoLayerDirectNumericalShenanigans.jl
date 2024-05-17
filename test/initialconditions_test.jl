# Need to write these to check things are set at the right level in the right field
diffusivities = (ν = 1e-4, κ = (S = 1e-5, T = 1e-5))
resolution = (Nx = 10, Ny = 10, Nz = 500)
model = DNSModel(architecture, DOMAIN_EXTENT, resolution, diffusivities)

## set initial conditions, currently there are four options available in this submodule
initial_conditions = StableTwoLayerInitialConditions(0, 0, 0, 0, 0, 0)
profile_function = HyperbolicTangent(INTERFACE_LOCATION, 50.0)
z = znodes(model.grid, Center(), Center(), Center())
depth_idx = findall(rand(z) .== z)[1]

tracer_profile_perturbations = (SalinityGaussianProfile(z[depth_idx], 0.0, 1.5),
                                TemperatureGaussianProfile(z[depth_idx], 0.0, 1.5))
tracer_blob_perturbations = (SalinityGaussianBlob(z[depth_idx], [0.0, 0.0], 1.5),
                             TemperatureGaussianBlob(z[depth_idx], [0.0, 0.0], 1.5))
tracer_noise_perturbations = (SalinityNoise(z[depth_idx], 1.0), TemperatureNoise(z[depth_idx], 1.0))
tracer_noise_perturbations_vec = (SalinityNoise(z[depth_idx-1:depth_idx+1], fill(1.0, 3)),
                                  TemperatureNoise(z[depth_idx-1:depth_idx+1], fill(1.0, 3)))
tracer_profile_noise_perturbations = ((SalinityGaussianProfile(z[depth_idx], 0.0, 1.5),
                                       SalinityNoise(z[depth_idx], 1.0)),
                                       (TemperatureGaussianProfile(z[depth_idx], 0.0, 1.5),
                                        TemperatureNoise(z[depth_idx], 1.0)))

function tracer_profile(dns::TwoLayerDNS)

    S, T = interior(dns.model.tracers.S, :, :, :), interior(dns.model.tracers.T, :, :, :)

    S_profile, T_profile = false, false
    if dns.tracer_perturbation isa SalinityGaussianProfile
        T_profile = isempty(findall(T .!= 0))
        S_profile = S[rand(1:resolution.Nx), rand(1:resolution.Ny), :] .==
                    [perturb_tracer(_z, dns.tracer_perturbation) for _z in z]
    elseif dns.tracer_perturbation isa TemperatureGaussianProfile
        T_profile = T[rand(1:resolution.Nx), rand(1:resolution.Ny), :] .==
                    [perturb_tracer(_z, dns.tracer_perturbation) for _z in z]
        S_profile = isempty(findall(S .!= 0))
    end

    return (T_profile, S_profile)

end
function tracer_blob(dns::TwoLayerDNS)

    S, T = interior(dns.model.tracers.S, :, :, :), interior(dns.model.tracers.T, :, :, :)

    find_T, find_S = false, false
    if dns.tracer_perturbation isa SalinityGaussianBlob
        find_T = isempty(findall(T .!= 0))
        find_S = findall(S .!= 0)[1][3] == depth_idx
    elseif dns.tracer_perturbation isa TemperatureGaussianBlob
        find_T = findall(T .!= 0)[1][3] == depth_idx
        find_S = isempty(findall(S .!= 0))
    end

    return (find_T, find_S)

end
function tracer_noise(dns::TwoLayerDNS)

    S, T = interior(dns.model.tracers.S, :, :, :), interior(dns.model.tracers.T, :, :, :)

    find_T, find_S = false, false
    if dns.initial_noise isa SalinityNoise{<:Number}
        find_T = isempty(findall(T .!= 0))
        find_S = findall(S .!= 0)[1][3] == depth_idx
    elseif dns.initial_noise isa TemperatureNoise{<:Number}
        find_T = findall(T .!= 0)[1][3] == depth_idx
        find_S = isempty(findall(S .!= 0))
    elseif dns.initial_noise isa SalinityNoise{<:Vector}
        find_T = isempty(findall(T .!= 0))
        find_vec = findall(S .!= 0)
        find_S = [fv[3] for fv ∈ find_vec[1:100:end]] == depth_idx-1:depth_idx+1
    elseif dns.initial_noise isa TemperatureNoise{<:Vector}
        find_vec = findall(T .!= 0)
        find_T = [fv[3] for fv ∈ find_vec[1:100:end]] == depth_idx-1:depth_idx+1
        find_S = isempty(findall(S .!= 0))
    end

    return (find_T, find_S)

end
function tracer_stepchange(dns::TwoLayerDNS)

    S, T = interior(dns.model.tracers.S, :, :, :), interior(dns.model.tracers.T, :, :, :)

    S_upper = unique(S[1, 1, 1:depth_idx])[1] == dns.initial_conditions.S₀ˡ
    S_lower = unique(S[1, 1, depth_idx + 1:length(z)])[1] == dns.initial_conditions.S₀ᵘ
    T_upper = unique(T[1, 1, 1:depth_idx])[1] == dns.initial_conditions.T₀ˡ
    T_lower = unique(T[1, 1, depth_idx + 1:length(z)])[1] == dns.initial_conditions.T₀ᵘ

    return S_upper, S_lower, T_upper, T_lower

end
function tracer_stepchangelineargradient(dns::TwoLayerDNS, z)

    il = dns.profile_function.interface_location
    il_idx = findfirst(z .> il)
    depth = z[il_idx ]
    depth2 = z[il_idx  - 1]
    dSdz, dTdz = -dns.profile_function.dSdz, dns.profile_function.dTdz

    S₀ˡ, S₀ᵘ = dns.initial_conditions.S₀ˡ, dns.initial_conditions.S₀ᵘ
    S₀ˡ_offset, S₀ᵘ_offset = S₀ˡ - dSdz * depth2, S₀ᵘ - dSdz * depth
    S_upper = S₀ᵘ_offset .+ dSdz .* z[z .> il]
    S_lower = S₀ˡ_offset .+ dSdz .* z[z .≤ il]
    S = vcat(S_lower, S_upper)

    T₀ˡ, T₀ᵘ = dns.initial_conditions.T₀ˡ, dns.initial_conditions.T₀ᵘ
    T₀ˡ_offset, T₀ᵘ_offset = T₀ˡ - dTdz * depth2, T₀ᵘ - dTdz * depth
    T_upper = T₀ᵘ_offset .+ dTdz .* z[z .> il]
    T_lower = T₀ˡ_offset .+ dTdz .* z[z .≤ il]
    T = vcat(T_lower, T_upper)

    return S, T

end
