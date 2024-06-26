"""
    function find_depth(model::Oceananigans.AbstractModel)
Find first instance of `depth` in the `z` direction of the model grid. If the exact depth
is not found the first instance that is larger (i.e. deeper) than `depth` is returned.
"""
function find_depth(model::Oceananigans.AbstractModel, depth::Number)

    found_depth = model.architecture isa CPU ? begin
                                                z = znodes(model.grid, Center(), Center(), Center())
                                                depth_idx = findfirst(z .≥ depth)
                                                z[depth_idx]
                                               end :
                                                allowscalar() do
                                                z = znodes(model.grid, Center(), Center(), Center())
                                                depth_idx = findfirst(z .≥ depth)
                                                z[depth_idx]
                                               end

    return found_depth
end
function find_depth(model::Oceananigans.AbstractModel, depth::AbstractVector)

    sort!(depth)
    found_depth = model.architecture isa CPU ? begin
                                                z = znodes(model.grid, Center(), Center(), Center())
                                                depth_idx = findall(depth[1] .≤ z .≤ depth[2])
                                                z[depth_idx]
                                               end :
                                                allowscalar() do
                                                z = znodes(model.grid, Center(), Center(), Center())
                                                depth_idx = findall(depth[1] .≤ z .≤ depth[2])
                                                Array(z[depth_idx])
                                               end

    return found_depth
end
"""
    function set_two_layer_initial_conditions(dns::TwoLayerDNS)

Set initial conditions for a `TwoLayerDNS` that are smooth according to the
`profile_funciton` with or without a `tracer_perturbation` in the upper layer.
"""
function set_two_layer_initial_conditions!(dns::TwoLayerDNS)

    model, profile_function, initial_conditions, tracer_perturbation, initial_noise = dns

    set_two_layer_initial_conditions!(model, initial_conditions, profile_function,
                                      tracer_perturbation, initial_noise)

    return nothing

end
####
#### `Erf` methods
####
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::Erf,
                                           tracer_perturbation::Nothing,
                                           initial_noise::Nothing)

    κₛ, κₜ = model.closure.κ.S, model.closure.κ.T
    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = erf_tracer_solution(z, S₀, ΔS, κₛ, profile_function)
    initial_T_profile(x, y, z) = erf_tracer_solution(z, T₀, ΔT, κₜ, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
####
#### Salinity perturbations
####
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::Erf,
                                           tracer_perturbation::SalinityGaussianProfile,
                                           initial_noise::Nothing)

    κₛ, κₜ = model.closure.κ.S, model.closure.κ.T
    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = erf_tracer_solution(z, S₀, ΔS, κₛ, profile_function) +
                                    perturb_tracer(z, tracer_perturbation)

    initial_T_profile(x, y, z) = erf_tracer_solution(z, T₀, ΔT, κₜ, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::Erf,
                                           tracer_perturbation::SalinityGaussianBlob,
                                           initial_noise::Nothing)

    κₛ, κₜ = model.closure.κ.S, model.closure.κ.T
    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = erf_tracer_solution(z, S₀, ΔS, κₛ, profile_function) +
                                    perturb_tracer(x, y, z, tracer_perturbation)

    initial_T_profile(x, y, z) = erf_tracer_solution(z, T₀, ΔT, κₜ, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
####
#### Salinity + noise
####
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::Erf,
                                           tracer_perturbation::SalinityGaussianProfile,
                                           initial_noise::SalinityNoise)

    κₛ, κₜ = model.closure.κ.S, model.closure.κ.T
    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = erf_tracer_solution(z, S₀, ΔS, κₛ, profile_function) +
                                 perturb_tracer(z, tracer_perturbation) +
                                 perturb_tracer(z, initial_noise)

    initial_T_profile(x, y, z) = erf_tracer_solution(z, T₀, ΔT, κₜ, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::Erf,
                                           tracer_perturbation::SalinityGaussianProfile,
                                           initial_noise::VelocityNoise)

    κₛ, κₜ = model.closure.κ.S, model.closure.κ.T
    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = erf_tracer_solution(z, S₀, ΔS, κₛ, profile_function) +
                                 perturb_tracer(z, tracer_perturbation)

    initial_T_profile(x, y, z) = erf_tracer_solution(z, T₀, ΔT, κₜ, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    perturb_velocity!(model, initial_noise)

    return nothing

end
####
#### Temperature perturbations
####
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::Erf,
                                           tracer_perturbation::TemperatureGaussianProfile,
                                           initial_noise::Nothing)

    κₛ, κₜ = model.closure.κ.S, model.closure.κ.T
    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = erf_tracer_solution(z, S₀, ΔS, κₛ, profile_function)

    initial_T_profile(x, y, z) = erf_tracer_solution(z, T₀, ΔT, κₜ, profile_function) +
                                 perturb_tracer(z, tracer_perturbation)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::Erf,
                                           tracer_perturbation::TemperatureGaussianBlob,
                                           initial_noise::Nothing)

    κₛ, κₜ = model.closure.κ.S, model.closure.κ.T
    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = erf_tracer_solution(z, S₀, ΔS, κₛ, profile_function)

    initial_T_profile(x, y, z) = erf_tracer_solution(z, T₀, ΔT, κₜ, profile_function) +
                                 perturb_tracer(x, y, z, tracer_perturbation)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
####
#### Temperature + noise
####
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::Erf,
                                           tracer_perturbation::TemperatureGaussianProfile,
                                           initial_noise::TemperatureNoise)

    κₛ, κₜ = model.closure.κ.S, model.closure.κ.T
    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = erf_tracer_solution(z, S₀, ΔS, κₛ, profile_function)

    initial_T_profile(x, y, z) = erf_tracer_solution(z, T₀, ΔT, κₜ, profile_function) +
                                 perturb_tracer(z, tracer_perturbation) +
                                 perturb_tracer(z, initial_noise)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::Erf,
                                           tracer_perturbation::TemperatureGaussianProfile,
                                           initial_noise::VelocityNoise)

    κₛ, κₜ = model.closure.κ.S, model.closure.κ.T
    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = erf_tracer_solution(z, S₀, ΔS, κₛ, profile_function)

    initial_T_profile(x, y, z) = erf_tracer_solution(z, T₀, ΔT, κₜ, profile_function) +
                                 perturb_tracer(z, tracer_perturbation)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    perturb_velocity!(model, initial_noise)

    return nothing

end
####
#### Noise only
####
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::Erf,
                                           tracer_perturbation::Nothing,
                                           initial_noise::SalinityNoise)

    κₛ, κₜ = model.closure.κ.S, model.closure.κ.T
    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = erf_tracer_solution(z, S₀, ΔS, κₛ, profile_function) +
                                 perturb_tracer(z, initial_noise)

    initial_T_profile(x, y, z) = erf_tracer_solution(z, T₀, ΔT, κₜ, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::Erf,
                                           tracer_perturbation::Nothing,
                                           initial_noise::TemperatureNoise)

    κₛ, κₜ = model.closure.κ.S, model.closure.κ.T
    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = erf_tracer_solution(z, S₀, ΔS, κₛ, profile_function)

    initial_T_profile(x, y, z) = erf_tracer_solution(z, T₀, ΔT, κₜ, profile_function) +
                                 perturb_tracer(z, initial_noise)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::Erf,
                                           tracer_perturbation::Nothing,
                                           initial_noise::VelocityNoise)

    κₛ, κₜ = model.closure.κ.S, model.closure.κ.T
    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = erf_tracer_solution(z, S₀, ΔS, κₛ, profile_function)
    initial_T_profile(x, y, z) = erf_tracer_solution(z, T₀, ΔT, κₜ, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    perturb_velocity!(model, initial_noise)

    return nothing

end
####
#### `HyperbolicTangent` methods
####
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::HyperbolicTangent,
                                           tracer_perturbation::Nothing,
                                           initial_noise::Nothing)

    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = tanh_initial_condition(z, S₀, ΔS, profile_function)
    initial_T_profile(x, y, z) = tanh_initial_condition(z, T₀, ΔT, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
####
#### Salinity
####
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::HyperbolicTangent,
                                           tracer_perturbation::SalinityGaussianProfile,
                                           initial_noise::Nothing)

    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = tanh_initial_condition(z, S₀, ΔS, profile_function) +
                                 perturb_tracer(z, tracer_perturbation)

    initial_T_profile(x, y, z) = tanh_initial_condition(z, T₀, ΔT, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::HyperbolicTangent,
                                           tracer_perturbation::SalinityGaussianBlob,
                                           initial_noise::Nothing)

    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = tanh_initial_condition(z, S₀, ΔS, profile_function) +
                                    perturb_tracer(x, y, z, tracer_perturbation)

    initial_T_profile(x, y, z) = tanh_initial_condition(z, T₀, ΔT, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
####
#### Salinity + noise
#####
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::HyperbolicTangent,
                                           tracer_perturbation::SalinityGaussianProfile,
                                           initial_noise::SalinityNoise)

    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = tanh_initial_condition(z, S₀, ΔS, profile_function) +
                                 perturb_tracer(z, tracer_perturbation) +
                                 perturb_tracer(z, initial_noise)

    initial_T_profile(x, y, z) = tanh_initial_condition(z, T₀, ΔT, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::HyperbolicTangent,
                                           tracer_perturbation::SalinityGaussianProfile,
                                           initial_noise::VelocityNoise)

    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = tanh_initial_condition(z, S₀, ΔS, profile_function) +
                                 perturb_tracer(z, tracer_perturbation)

    initial_T_profile(x, y, z) = tanh_initial_condition(z, T₀, ΔT, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    perturb_velocity!(model, initial_noise)

    return nothing

end
####
#### Temperature
####
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::HyperbolicTangent,
                                           tracer_perturbation::TemperatureGaussianProfile,
                                           initial_noise::Nothing)

    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = tanh_initial_condition(z, S₀, ΔS, profile_function)

    initial_T_profile(x, y, z) = tanh_initial_condition(z, T₀, ΔT, profile_function) +
                                 perturb_tracer(z, tracer_perturbation)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::HyperbolicTangent,
                                           tracer_perturbation::TemperatureGaussianBlob,
                                           initial_noise::Nothing)

    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = tanh_initial_condition(z, S₀, ΔS, profile_function)

    initial_T_profile(x, y, z) = tanh_initial_condition(z, T₀, ΔT, profile_function) +
                                 perturb_tracer(x, y, z, tracer_perturbation)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
####
#### Temperature + noise
####
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::HyperbolicTangent,
                                           tracer_perturbation::TemperatureGaussianProfile,
                                           initial_noise::TemperatureNoise)

    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = tanh_initial_condition(z, S₀, ΔS, profile_function)

    initial_T_profile(x, y, z) = tanh_initial_condition(z, T₀, ΔT, profile_function) +
                                 perturb_tracer(z, tracer_perturbation) +
                                 perturb_tracer(z, initial_noise)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::HyperbolicTangent,
                                           tracer_perturbation::TemperatureGaussianProfile,
                                           initial_noise::VelocityNoise)

    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = tanh_initial_condition(z, S₀, ΔS, profile_function)

    initial_T_profile(x, y, z) = tanh_initial_condition(z, T₀, ΔT, profile_function) +
                                 perturb_tracer(z, tracer_perturbation)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    perturb_velocity!(model, initial_noise)

    return nothing

end
####
#### Noise
####
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::HyperbolicTangent,
                                           tracer_perturbation::Nothing,
                                           initial_noise::SalinityNoise)

    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = tanh_initial_condition(z, S₀, ΔS, profile_function) +
                                    perturb_tracer(z, initial_noise)

    initial_T_profile(x, y, z) = tanh_initial_condition(z, T₀, ΔT, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::HyperbolicTangent,
                                           tracer_perturbation::Nothing,
                                           initial_noise::TemperatureNoise)

    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = tanh_initial_condition(z, S₀, ΔS, profile_function)

    initial_T_profile(x, y, z) = tanh_initial_condition(z, T₀, ΔT, profile_function) +
                                 perturb_tracer(z, initial_noise)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                          initial_conditions::TwoLayerInitialConditions,
                                          profile_function::HyperbolicTangent,
                                          tracer_perturbation::Nothing,
                                          initial_noise::VelocityNoise)

    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = tanh_initial_condition(z, S₀, ΔS, profile_function)
    initial_T_profile(x, y, z) = tanh_initial_condition(z, T₀, ΔT, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    perturb_velocity!(model, initial_noise)

    return nothing

end
####
#### `MidPoint` methods
####
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::MidPoint,
                                           tracer_perturbation::Nothing,
                                           initial_noise::Nothing)

    S₀ᵘ, S₀ˡ = initial_conditions.S₀ᵘ, initial_conditions.S₀ˡ
    T₀ᵘ, T₀ˡ = initial_conditions.T₀ᵘ, initial_conditions.T₀ˡ

    initial_S_profile(x, y, z) = midpoint(z, S₀ˡ, S₀ᵘ, profile_function)
    initial_T_profile(x, y, z) = midpoint(z, T₀ˡ, T₀ᵘ, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
####
#### Salinity
####
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::MidPoint,
                                           tracer_perturbation::SalinityGaussianProfile,
                                           initial_noise::Nothing)

    S₀ᵘ, S₀ˡ = initial_conditions.S₀ᵘ, initial_conditions.S₀ˡ
    T₀ᵘ, T₀ˡ = initial_conditions.T₀ᵘ, initial_conditions.T₀ˡ

    initial_S_profile(x, y, z) = midpoint(z, S₀ˡ, S₀ᵘ, profile_function) +
                                 perturb_tracer(z, tracer_perturbation)
    initial_T_profile(x, y, z) = midpoint(z, T₀ˡ, T₀ᵘ, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::MidPoint,
                                           tracer_perturbation::SalinityGaussianBlob,
                                           initial_noise::Nothing)

    S₀ᵘ, S₀ˡ = initial_conditions.S₀ᵘ, initial_conditions.S₀ˡ
    T₀ᵘ, T₀ˡ = initial_conditions.T₀ᵘ, initial_conditions.T₀ˡ

    initial_S_profile(x, y, z) = midpoint(z, S₀ˡ, S₀ᵘ, profile_function) +
                                 perturb_tracer(x, y, z, tracer_perturbation)
    initial_T_profile(x, y, z) = midpoint(z, T₀ˡ, T₀ᵘ, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
####
#### Salinity + noise
#####
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::MidPoint,
                                           tracer_perturbation::SalinityGaussianProfile,
                                           initial_noise::SalinityNoise)

    S₀ᵘ, S₀ˡ = initial_conditions.S₀ᵘ, initial_conditions.S₀ˡ
    T₀ᵘ, T₀ˡ = initial_conditions.T₀ᵘ, initial_conditions.T₀ˡ

    initial_S_profile(x, y, z) = midpoint(z, S₀ˡ, S₀ᵘ, profile_function) +
                                 perturb_tracer(z, tracer_perturbation) +
                                 perturb_tracer(z, initial_noise)
    initial_T_profile(x, y, z) = midpoint(z, T₀ˡ, T₀ᵘ, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::MidPoint,
                                           tracer_perturbation::SalinityGaussianProfile,
                                           initial_noise::VelocityNoise)

    S₀ᵘ, S₀ˡ = initial_conditions.S₀ᵘ, initial_conditions.S₀ˡ
    T₀ᵘ, T₀ˡ = initial_conditions.T₀ᵘ, initial_conditions.T₀ˡ

    initial_S_profile(x, y, z) = midpoint(z, S₀ˡ, S₀ᵘ, profile_function) +
                                 perturb_tracer(z, tracer_perturbation)
    initial_T_profile(x, y, z) = midpoint(z, T₀ˡ, T₀ᵘ, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    perturb_velocity!(model, initial_noise)

    return nothing

end
####
#### Temperature
####
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::MidPoint,
                                           tracer_perturbation::TemperatureGaussianProfile,
                                           initial_noise::Nothing)

    S₀ᵘ, S₀ˡ = initial_conditions.S₀ᵘ, initial_conditions.S₀ˡ
    T₀ᵘ, T₀ˡ = initial_conditions.T₀ᵘ, initial_conditions.T₀ˡ

    initial_S_profile(x, y, z) = midpoint(z, S₀ˡ, S₀ᵘ, profile_function)
    initial_T_profile(x, y, z) = midpoint(z, T₀ˡ, T₀ᵘ, profile_function) +
                                 perturb_tracer(z, tracer_perturbation)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::MidPoint,
                                           tracer_perturbation::TemperatureGaussianBlob,
                                           initial_noise::Nothing)

    S₀ᵘ, S₀ˡ = initial_conditions.S₀ᵘ, initial_conditions.S₀ˡ
    T₀ᵘ, T₀ˡ = initial_conditions.T₀ᵘ, initial_conditions.T₀ˡ

    initial_S_profile(x, y, z) = midpoint(z, S₀ˡ, S₀ᵘ, profile_function)
    initial_T_profile(x, y, z) = midpoint(z, T₀ˡ, T₀ᵘ, profile_function) +
                                 perturb_tracer(x, y, z, tracer_perturbation)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
####
#### Temperature + noise
####
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::MidPoint,
                                           tracer_perturbation::TemperatureGaussianProfile,
                                           initial_noise::TemperatureNoise)

    S₀ᵘ, S₀ˡ = initial_conditions.S₀ᵘ, initial_conditions.S₀ˡ
    T₀ᵘ, T₀ˡ = initial_conditions.T₀ᵘ, initial_conditions.T₀ˡ

    initial_S_profile(x, y, z) = midpoint(z, S₀ˡ, S₀ᵘ, profile_function)
    initial_T_profile(x, y, z) = midpoint(z, T₀ˡ, T₀ᵘ, profile_function) +
                                 perturb_tracer(z, tracer_perturbation) +
                                 perturb_tracer(z, initial_noise)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::MidPoint,
                                           tracer_perturbation::TemperatureGaussianProfile,
                                           initial_noise::VelocityNoise)

    S₀ᵘ, S₀ˡ = initial_conditions.S₀ᵘ, initial_conditions.S₀ˡ
    T₀ᵘ, T₀ˡ = initial_conditions.T₀ᵘ, initial_conditions.T₀ˡ

    initial_S_profile(x, y, z) = midpoint(z, S₀ˡ, S₀ᵘ, profile_function)
    initial_T_profile(x, y, z) = midpoint(z, T₀ˡ, T₀ᵘ, profile_function) +
                                perturb_tracer(z, tracer_perturbation)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    perturb_velocity!(model, initial_noise)

    return nothing

end
####
#### Noise
####
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::MidPoint,
                                           tracer_perturbation::Nothing,
                                           initial_noise::SalinityNoise)

    S₀ᵘ, S₀ˡ = initial_conditions.S₀ᵘ, initial_conditions.S₀ˡ
    T₀ᵘ, T₀ˡ = initial_conditions.T₀ᵘ, initial_conditions.T₀ˡ

    initial_S_profile(x, y, z) = midpoint(z, S₀ˡ, S₀ᵘ, profile_function) +
                                 perturb_tracer(z, initial_noise)
    initial_T_profile(x, y, z) = midpoint(z, T₀ˡ, T₀ᵘ, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::MidPoint,
                                           tracer_perturbation::Nothing,
                                           initial_noise::TemperatureNoise)

    S₀ᵘ, S₀ˡ = initial_conditions.S₀ᵘ, initial_conditions.S₀ˡ
    T₀ᵘ, T₀ˡ = initial_conditions.T₀ᵘ, initial_conditions.T₀ˡ

    initial_S_profile(x, y, z) = midpoint(z, S₀ˡ, S₀ᵘ, profile_function)
    initial_T_profile(x, y, z) = midpoint(z, T₀ˡ, T₀ᵘ, profile_function) +
                                 perturb_tracer(z, initial_noise)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::MidPoint,
                                           tracer_perturbation::Nothing,
                                           initial_noise::VelocityNoise)

    S₀ᵘ, S₀ˡ = initial_conditions.S₀ᵘ, initial_conditions.S₀ˡ
    T₀ᵘ, T₀ˡ = initial_conditions.T₀ᵘ, initial_conditions.T₀ˡ

    initial_S_profile(x, y, z) = midpoint(z, S₀ˡ, S₀ᵘ, profile_function)
    initial_T_profile(x, y, z) = midpoint(z, T₀ˡ, T₀ᵘ, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

    perturb_velocity!(model, initial_noise)

end
####
#### Step Change methods
####
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::StepChange,
                                           tracer_perturbation::Nothing,
                                           initial_noise::Nothing)

    S₀ᵘ, S₀ˡ = initial_conditions.S₀ᵘ, initial_conditions.S₀ˡ
    T₀ᵘ, T₀ˡ = initial_conditions.T₀ᵘ, initial_conditions.T₀ˡ

    initial_S_profile(x, y, z) = Heaviside(z, S₀ˡ, S₀ᵘ, profile_function)
    initial_T_profile(x, y, z) = Heaviside(z, T₀ˡ, T₀ᵘ, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
####
#### Salinity
####
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::StepChange,
                                           tracer_perturbation::SalinityGaussianProfile,
                                           initial_noise::Nothing)

    S₀ᵘ, S₀ˡ = initial_conditions.S₀ᵘ, initial_conditions.S₀ˡ
    T₀ᵘ, T₀ˡ = initial_conditions.T₀ᵘ, initial_conditions.T₀ˡ

    initial_S_profile(x, y, z) = Heaviside(z, S₀ˡ, S₀ᵘ, profile_function) +
                                 perturb_tracer(z, tracer_perturbation)
    initial_T_profile(x, y, z) = Heaviside(z, T₀ˡ, T₀ᵘ, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::StepChange,
                                           tracer_perturbation::SalinityGaussianBlob,
                                           initial_noise::Nothing)

    S₀ᵘ, S₀ˡ = initial_conditions.S₀ᵘ, initial_conditions.S₀ˡ
    T₀ᵘ, T₀ˡ = initial_conditions.T₀ᵘ, initial_conditions.T₀ˡ

    initial_S_profile(x, y, z) = Heaviside(z, S₀ˡ, S₀ᵘ, profile_function) +
                                 perturb_tracer(x, y, z, tracer_perturbation)
    initial_T_profile(x, y, z) = Heaviside(z, T₀ˡ, T₀ᵘ, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
####
#### Salinity + noise
#####
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::StepChange,
                                           tracer_perturbation::SalinityGaussianProfile,
                                           initial_noise::SalinityNoise)

    S₀ᵘ, S₀ˡ = initial_conditions.S₀ᵘ, initial_conditions.S₀ˡ
    T₀ᵘ, T₀ˡ = initial_conditions.T₀ᵘ, initial_conditions.T₀ˡ

    initial_S_profile(x, y, z) = Heaviside(z, S₀ˡ, S₀ᵘ, profile_function) +
                                 perturb_tracer(z, tracer_perturbation) +
                                 perturb_tracer(z, initial_noise)
    initial_T_profile(x, y, z) = Heaviside(z, T₀ˡ, T₀ᵘ, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::StepChange,
                                           tracer_perturbation::SalinityGaussianProfile,
                                           initial_noise::VelocityNoise)

    S₀ᵘ, S₀ˡ = initial_conditions.S₀ᵘ, initial_conditions.S₀ˡ
    T₀ᵘ, T₀ˡ = initial_conditions.T₀ᵘ, initial_conditions.T₀ˡ

    initial_S_profile(x, y, z) = Heaviside(z, S₀ˡ, S₀ᵘ, profile_function) +
                                 perturb_tracer(z, tracer_perturbation)
    initial_T_profile(x, y, z) = Heaviside(z, T₀ˡ, T₀ᵘ, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    perturb_velocity!(model, initial_noise)

    return nothing

end
####
#### Temperature
####
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::StepChange,
                                           tracer_perturbation::TemperatureGaussianProfile,
                                           initial_noise::Nothing)

    S₀ᵘ, S₀ˡ = initial_conditions.S₀ᵘ, initial_conditions.S₀ˡ
    T₀ᵘ, T₀ˡ = initial_conditions.T₀ᵘ, initial_conditions.T₀ˡ

    initial_S_profile(x, y, z) = Heaviside(z, S₀ˡ, S₀ᵘ, profile_function)
    initial_T_profile(x, y, z) = Heaviside(z, T₀ˡ, T₀ᵘ, profile_function) +
                                 perturb_tracer(z, tracer_perturbation)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::StepChange,
                                           tracer_perturbation::TemperatureGaussianBlob,
                                           initial_noise::Nothing)

    S₀ᵘ, S₀ˡ = initial_conditions.S₀ᵘ, initial_conditions.S₀ˡ
    T₀ᵘ, T₀ˡ = initial_conditions.T₀ᵘ, initial_conditions.T₀ˡ

    initial_S_profile(x, y, z) = Heaviside(z, S₀ˡ, S₀ᵘ, profile_function)
    initial_T_profile(x, y, z) = Heaviside(z, T₀ˡ, T₀ᵘ, profile_function) +
                                 perturb_tracer(x, y, z, tracer_perturbation)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
####
#### Temperature + noise
####
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::StepChange,
                                           tracer_perturbation::TemperatureGaussianProfile,
                                           initial_noise::TemperatureNoise)

    S₀ᵘ, S₀ˡ = initial_conditions.S₀ᵘ, initial_conditions.S₀ˡ
    T₀ᵘ, T₀ˡ = initial_conditions.T₀ᵘ, initial_conditions.T₀ˡ

    initial_S_profile(x, y, z) = Heaviside(z, S₀ˡ, S₀ᵘ, profile_function)
    initial_T_profile(x, y, z) = Heaviside(z, T₀ˡ, T₀ᵘ, profile_function) +
                                 perturb_tracer(z, tracer_perturbation) +
                                 perturb_tracer(z, initial_noise)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::StepChange,
                                           tracer_perturbation::TemperatureGaussianProfile,
                                           initial_noise::VelocityNoise)

    S₀ᵘ, S₀ˡ = initial_conditions.S₀ᵘ, initial_conditions.S₀ˡ
    T₀ᵘ, T₀ˡ = initial_conditions.T₀ᵘ, initial_conditions.T₀ˡ

    initial_S_profile(x, y, z) = Heaviside(z, S₀ˡ, S₀ᵘ, profile_function)
    initial_T_profile(x, y, z) = Heaviside(z, T₀ˡ, T₀ᵘ, profile_function) +
                                 perturb_tracer(z, tracer_perturbation)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    perturb_velocity!(model, initial_noise)

    return nothing

end
####
#### Noise
####
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::StepChange,
                                           tracer_perturbation::Nothing,
                                           initial_noise::SalinityNoise)

    S₀ᵘ, S₀ˡ = initial_conditions.S₀ᵘ, initial_conditions.S₀ˡ
    T₀ᵘ, T₀ˡ = initial_conditions.T₀ᵘ, initial_conditions.T₀ˡ

    initial_S_profile(x, y, z) = Heaviside(z, S₀ˡ, S₀ᵘ, profile_function) +
                                 perturb_tracer(z, initial_noise)
    initial_T_profile(x, y, z) = Heaviside(z, T₀ˡ, T₀ᵘ, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::StepChange,
                                           tracer_perturbation::Nothing,
                                           initial_noise::TemperatureNoise)

    S₀ᵘ, S₀ˡ = initial_conditions.S₀ᵘ, initial_conditions.S₀ˡ
    T₀ᵘ, T₀ˡ = initial_conditions.T₀ᵘ, initial_conditions.T₀ˡ

    initial_S_profile(x, y, z) = Heaviside(z, S₀ˡ, S₀ᵘ, profile_function)
    initial_T_profile(x, y, z) = Heaviside(z, T₀ˡ, T₀ᵘ, profile_function) +
                                 perturb_tracer(z, initial_noise)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::StepChange,
                                           tracer_perturbation::Nothing,
                                           initial_noise::VelocityNoise)

    S₀ᵘ, S₀ˡ = initial_conditions.S₀ᵘ, initial_conditions.S₀ˡ
    T₀ᵘ, T₀ˡ = initial_conditions.T₀ᵘ, initial_conditions.T₀ˡ

    initial_S_profile(x, y, z) = Heaviside(z, S₀ˡ, S₀ᵘ, profile_function)
    initial_T_profile(x, y, z) = Heaviside(z, T₀ˡ, T₀ᵘ, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

    perturb_velocity!(model, initial_noise)

end
####
#### Step Change with linear gradient
####
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                            initial_conditions::TwoLayerInitialConditions,
                                            profile_function::StepChangeLinearGradient,
                                            tracer_perturbation::Nothing,
                                            initial_noise::Nothing)

    S₀ᵘ, S₀ˡ = initial_conditions.S₀ᵘ, initial_conditions.S₀ˡ
    T₀ᵘ, T₀ˡ = initial_conditions.T₀ᵘ, initial_conditions.T₀ˡ

    initial_S_profile(x, y, z) = Heaviside_with_linear_gradient(z, S₀ˡ, S₀ᵘ, profile_function)
    initial_T_profile(x, y, z) = Heaviside_with_linear_gradient(z, T₀ˡ, T₀ᵘ, profile_function, tracer = :T)

    set!(model, S = initial_S_profile, T = initial_T_profile)

return nothing

end
#### Step Change with linear gradient + salinity noise
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                            initial_conditions::TwoLayerInitialConditions,
                                            profile_function::StepChangeLinearGradient,
                                            tracer_perturbation::Nothing,
                                            initial_noise::SalinityNoise)

S₀ᵘ, S₀ˡ = initial_conditions.S₀ᵘ, initial_conditions.S₀ˡ
T₀ᵘ, T₀ˡ = initial_conditions.T₀ᵘ, initial_conditions.T₀ˡ

initial_S_profile(x, y, z) = Heaviside_with_linear_gradient(z, S₀ˡ, S₀ᵘ, profile_function) +
                             perturb_tracer(z, initial_noise)
initial_T_profile(x, y, z) = Heaviside_with_linear_gradient(z, T₀ˡ, T₀ᵘ, profile_function, tracer = :T)

set!(model, S = initial_S_profile, T = initial_T_profile)

return nothing

end
