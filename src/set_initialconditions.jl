"""
    function set_two_layer_initial_conditions(dns::TwoLayerDNS)

Set initial conditions for a `TwoLayerDNS` that are smooth according to the
`profile_funciton` with or without a `salinity_perturbation` in the upper layer.
"""
function set_two_layer_initial_conditions!(dns::TwoLayerDNS)

    model, initial_conditions, profile_function, salinity_perturbation =
        dns.model, dns.initial_conditions, dns.profile_function, dns.salinity_perturbation

    isnothing(salinity_perturbation) ? set_two_layer_initial_conditions!(model,
                                                                         initial_conditions,
                                                                         profile_function) :
                                       set_two_layer_initial_conditions!(model,
                                                                         initial_conditions,
                                                                         profile_function,
                                                                         salinity_perturbation)
    return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::Erf)

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
                                           salinity_perturbation::GaussianProfile)

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
                                           salinity_noise::RandomPerturbations)

κₛ, κₜ = model.closure.κ.S, model.closure.κ.T
S₀ = initial_conditions.S₀ˡ
ΔS = initial_conditions.ΔS₀
T₀ = initial_conditions.T₀ˡ
ΔT = initial_conditions.ΔT₀

initial_S_profile(x, y, z) = erf_tracer_solution(z, S₀, ΔS, κₛ, profile_function) +
                                perturb_salinity(z, salinity_perturbation) +
                                perturb_salinity(z, salinity_noise)

initial_T_profile(x, y, z) = erf_tracer_solution(z, T₀, ΔT, κₜ, profile_function)

set!(model, S = initial_S_profile, T = initial_T_profile)

return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::Erf,
                                           salinity_perturbation::GaussianBlob)

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
                                           salinity_perturbation::RandomPerturbations)

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
                                           profile_function::HyperbolicTangent)

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
                                           salinity_perturbation::GaussianProfile)

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
                                           salinity_noise::RandomPerturbations)

S₀ = initial_conditions.S₀ˡ
ΔS = initial_conditions.ΔS₀
T₀ = initial_conditions.T₀ˡ
ΔT = initial_conditions.ΔT₀

initial_S_profile(x, y, z) = tanh_initial_condition(z, S₀, ΔS, profile_function) +
                                perturb_salinity(z, salinity_perturbation) +
                                perturb_salinity(z, salinity_noise)

initial_T_profile(x, y, z) = tanh_initial_condition(z, T₀, ΔT, profile_function)

set!(model, S = initial_S_profile, T = initial_T_profile)

return nothing

end
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           initial_conditions::TwoLayerInitialConditions,
                                           profile_function::HyperbolicTangent,
                                           salinity_perturbation::GaussianBlob)

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
                                           salinity_perturbation::RandomPerturbations)

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
"""
    function add_velocity_random_noise!(model::Oceananigans.AbstractModel, noise_magnitude::Number,
                                        location::Number; horizontal = false)
Add standard normally distributed random noise, scaled by `noise_magnitude`, to the
horizontal velocity fields at the interface of the upper and lower layers. If
`location` for where noise should be seeded is not provided then the noise is added
everywhere in the domain. To only add horizontal random noise (i.e. in the `u` and `v`
velocity fields) set `true` for the `horizontal` keyword argument.
"""
function add_velocity_random_noise!(model::Oceananigans.AbstractModel,
                                    noise_magnitude::Number, location::Number;
                                    horizontal = false)

    add_noise(x, y, z) = round(z; digits = 4) == location ? noise_magnitude * randn() : 0

    horizontal == true ? set!(model, u = add_noise, v = add_noise) :
                         set!(model, u = add_noise, v = add_noise, w = add_noise)

    return nothing

end
function add_velocity_random_noise!(model::Oceananigans.AbstractModel,
                                    noise_magnitude::Number; horizontal = false)

    add_noise(x, y, z) = noise_magnitude * randn()

    horizontal == true ? set!(model, u = add_noise, v = add_noise) :
                         set!(model, u = add_noise, v = add_noise, w = add_noise)

    return nothing

end
