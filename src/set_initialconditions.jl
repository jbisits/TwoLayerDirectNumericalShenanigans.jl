"""
    function set_two_layer_initial_conditions(dns::TwoLayerDNS)

Set initial conditions for a `TwoLayerDNS` that are smooth according to the
`profile_funciton` with or without a `salinity_perturbation` in the upper layer.
"""
function set_two_layer_initial_conditions!(dns::TwoLayerDNS)

    model, initial_conditions, profile_function, salinity_perturbation, initial_noise =
        dns.model, dns.initial_conditions, dns.profile_function, dns.salinity_perturbation, dns.initial_noise

    set_two_layer_initial_conditions!(model, initial_conditions, profile_function,
                                      salinity_perturbation, initial_noise)

    return nothing

end
###### `Erf` methods
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::Erf,
                                           salinity_perturbation::Nothing,
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
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::Erf,
                                           salinity_perturbation::GaussianProfile,
                                           initial_noise::Nothing)

    κₛ, κₜ = model.closure.κ.S, model.closure.κ.T
    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = erf_tracer_solution(z, S₀, ΔS, κₛ, profile_function) +
                                    perturb_salinity(z, salinity_perturbation)

    initial_T_profile(x, y, z) = erf_tracer_solution(z, T₀, ΔT, κₜ, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::Erf,
                                           salinity_perturbation::GaussianProfile,
                                           initial_noise::SalinityNoise)

    κₛ, κₜ = model.closure.κ.S, model.closure.κ.T
    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = erf_tracer_solution(z, S₀, ΔS, κₛ, profile_function) +
                                    perturb_salinity(z, salinity_perturbation) +
                                    perturb_salinity(z, initial_noise)

    initial_T_profile(x, y, z) = erf_tracer_solution(z, T₀, ΔT, κₜ, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::Erf,
                                           salinity_perturbation::GaussianProfile,
                                           initial_noise::VelocityNoise)

    κₛ, κₜ = model.closure.κ.S, model.closure.κ.T
    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = erf_tracer_solution(z, S₀, ΔS, κₛ, profile_function) +
                                 perturb_salinity(z, salinity_perturbation)

    initial_T_profile(x, y, z) = erf_tracer_solution(z, T₀, ΔT, κₜ, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    perturb_velocity!(model, initial_noise)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::Erf,
                                           salinity_perturbation::GaussianBlob,
                                           initial_noise::Nothing)

    κₛ, κₜ = model.closure.κ.S, model.closure.κ.T
    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = erf_tracer_solution(z, S₀, ΔS, κₛ, profile_function) +
                                    perturb_salinity(x, y, z, salinity_perturbation)

    initial_T_profile(x, y, z) = erf_tracer_solution(z, T₀, ΔT, κₜ, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::Erf,
                                           salinity_perturbation::Nothing,
                                           initial_noise::SalinityNoise)

    κₛ, κₜ = model.closure.κ.S, model.closure.κ.T
    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = erf_tracer_solution(z, S₀, ΔS, κₛ, profile_function) +
                                    perturb_salinity(z, initial_noise)

    initial_T_profile(x, y, z) = erf_tracer_solution(z, T₀, ΔT, κₜ, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
    initial_conditions::TwoLayerInitialConditions,
    profile_function::Erf,
    salinity_perturbation::Nothing,
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
###### `HyperbolicTangent` methods
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::HyperbolicTangent,
                                           salinity_perturbation::Nothing,
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
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::HyperbolicTangent,
                                           salinity_perturbation::GaussianProfile,
                                           initial_noise::Nothing)

    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = tanh_initial_condition(z, S₀, ΔS, profile_function) +
                                    perturb_salinity(z, salinity_perturbation)

    initial_T_profile(x, y, z) = tanh_initial_condition(z, T₀, ΔT, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::HyperbolicTangent,
                                           salinity_perturbation::GaussianProfile,
                                           initial_noise::SalinityNoise)

    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = tanh_initial_condition(z, S₀, ΔS, profile_function) +
                                    perturb_salinity(z, salinity_perturbation) +
                                    perturb_salinity(z, initial_noise)

    initial_T_profile(x, y, z) = tanh_initial_condition(z, T₀, ΔT, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::HyperbolicTangent,
                                           salinity_perturbation::GaussianProfile,
                                           initial_noise::VelocityNoise)

    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = tanh_initial_condition(z, S₀, ΔS, profile_function) +
                                 perturb_salinity(z, salinity_perturbation)

    initial_T_profile(x, y, z) = tanh_initial_condition(z, T₀, ΔT, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    perturb_velocity!(model, initial_noise)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::HyperbolicTangent,
                                           salinity_perturbation::GaussianBlob,
                                           initial_noise::Nothing)

    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = tanh_initial_condition(z, S₀, ΔS, profile_function) +
                                    perturb_salinity(x, y, z, salinity_perturbation)

    initial_T_profile(x, y, z) = tanh_initial_condition(z, T₀, ΔT, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::HyperbolicTangent,
                                           salinity_perturbation::Nothing,
                                           initial_noise::SalinityNoise)

    S₀ = initial_conditions.S₀ˡ
    ΔS = initial_conditions.ΔS₀
    T₀ = initial_conditions.T₀ˡ
    ΔT = initial_conditions.ΔT₀

    initial_S_profile(x, y, z) = tanh_initial_condition(z, S₀, ΔS, profile_function) +
                                    perturb_salinity(z, initial_noise)

    initial_T_profile(x, y, z) = tanh_initial_condition(z, T₀, ΔT, profile_function)

    set!(model, S = initial_S_profile, T = initial_T_profile)

    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                          initial_conditions::TwoLayerInitialConditions,
                                          profile_function::HyperbolicTangent,
                                          salinity_perturbation::Nothing,
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
