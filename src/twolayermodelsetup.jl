"""
    function set_two_layer_initial_conditions(model::AbstractModel, S::NamedTuple, Θ::NamedTuple)
Set initial conditions for temperature and salinity in a two layer model. Initial values for
the absolute salinity (`S`) and conservative temperature (`Θ`) in each layer must be passed
to the function as a `NamedTuple` in the form `S = (upper = value, lower = value)` as well
as the model where the initial conditions should be set. **Note** by default the interface
is in the middle of the domain, `z_centre = 0.5`, with width of `interface_thickness = 100`
These can be altered by passing keyworkd arguments.
"""
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel, S::NamedTuple,
                                           Θ::NamedTuple; z_centre = 0.5,
                                           interface_thickness = 100)

    ΔS = (S.upper - S.lower) / 2
    ΔΘ = (Θ.upper - Θ.lower) / 2

    initial_S_profile(x, y, z) = ΔS * tanh(interface_thickness * (z + z_centre)) + (S.lower + ΔS) #+ perturb_salintiy(z)
    initial_Θ_profile(x, y, z) = ΔΘ * tanh(interface_thickness * (z + z_centre)) + (Θ.lower + ΔΘ)

    set!(model, S = initial_S_profile, T = initial_Θ_profile)

    return nothing

end

"""
    function perturb_salintiy(z; salinity_pertubration)
Where and what value to add to perturb the salinity initial condition.
"""
function perturb_salintiy(z; salinity_pertubration = 0.032)
    if -0.2 > z > -0.25
        salinity_pertubration
    else
        0
    end
end

"""
    function viscosity_and_diffusivity(Prandtl_number::Number, Raleigh_number::Number)
Return the kinematic viscosity (`ν`) and the diffusivity (`κ`) from the
[Prandtl number](https://en.wikipedia.org/wiki/Prandtl_number) and
the [Raleigh number](https://en.wikipedia.org/wiki/Rayleigh_number).
The Prandtl number is
```math
    Pr = \\frac{ν}{κ}
```
and the Raleigh number is
```math
    Ra = \\frac{ΔρgL³}{μκ} = \\frac{bL³}{νκ}.
```
Using these expressions we can derive `ν` and `κ` in terms of the non-dimensional numbers,
the buoyancy and the domain size
```math
ν = \\sqrt{\\frac{PrbL³}{Ra}} \\quad κ = \\sqrt{\\frac{bL³}{RaPr}}.
```
"""
function viscosity_and_diffusivity(model::Oceananigans.AbstractModel,
                                   Prandtl_number::Number, Raleigh_number::Number)

    T, S = model.tracers.T, model.tracers.S

    #b = model.buoyancy

    return nothing
end

"""
    function simulation_progress(sim)
Useful progress messaging for simulation runs
"""
simulation_progress(sim) = @printf("i: % 6d, sim time: % 1.3f, wall time: % 10s, Δt: % 1.4f, advective CFL: %.2e, diffusive CFL: %.2e\n",
                                    iteration(sim), time(sim), prettytime(sim.run_wall_time),
                                    sim.Δt, AdvectiveCFL(sim.Δt)(sim.model),
                                    DiffusiveCFL(sim.Δt)(sim.model))
