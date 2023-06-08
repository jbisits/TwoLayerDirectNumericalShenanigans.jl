"""
    function DNS_cabbeling(resolution::Tuple, diffusivities::NamedTuple)
Setup a Direct Numerical Simulation in a 0.05m × 0.05m × 1m box on `architecture`
(`CPU()` or `GPU()`) with grid strecthing over approximately the bottom 3rd of the domain.
The argument `resoltuion` is a `NamedTuple` and needs to be passed as (Nx = , Ny = , Nz = )
and `diffusivities` is a also a `NamedTuple` of form `(ν = , κ = )` for the momentum and
tracer diffusivities respectively (salt and heat can have different diffusivities).
The equation of state used to evolve the buoyancy is the [polynomial approximation to the
full non-linear TEOS10 equation of state appropriate for Boussinesq models](https://clima.github.io/SeawaterPolynomials.jl/dev/API/#SeawaterPolynomials.TEOS10).
The keyword argument `reference_density` can be passed to evolve the density
about. If nothing is passed default reference density of 1020kgm⁻³ is used. The grid spacing
at the top of the domain is controlled by the keyword argument `refinement` while the rate
the bottom is stretched at is controlled by `stretching`.
"""
function DNS_cabbeling(architecture, resolution::NamedTuple, diffusivities::NamedTuple;
                       reference_density = nothing,
                       refinement = 1.2,
                       stretching = 100)

    Lx, Ly, Lz = 0.1, 0.1, 1
    Nx, Ny, Nz = resolution.Nx, resolution.Ny, resolution.Nz

    grid = RectilinearGrid(architecture,
                           topology = (Periodic, Periodic, Bounded),
                           size = (Nx, Ny, Nz),
                           x = (-Lx/2, Lx/2),
                           y = (-Ly/2, Ly/2),
                           z = grid_stretching(Lz, Nz, refinement, stretching))

    eos = isnothing(reference_density) ? SeawaterPolynomials.TEOS10EquationOfState() :
                            SeawaterPolynomials.TEOS10EquationOfState(; reference_density)
    buoyancy = SeawaterBuoyancy(equation_of_state = eos)
    tracers = (:T, :S)

    closure = ScalarDiffusivity(; ν = diffusivities.ν, κ = diffusivities.κ)

    timestepper = :RungeKutta3

    advection = CenteredSecondOrder()

    return NonhydrostaticModel(; grid, buoyancy, tracers, closure, timestepper, advection)

end
"""
    function grid_stretching(Lz, Nz; refinement, stretching)
Stretch the vertical coordinate of the grid. This function was taken from an
[Oceananigans.jl example](https://clima.github.io/OceananigansDocumentation/dev/generated/ocean_wind_mixing_and_convection/).
The rate of stretching at the bottom is controlled by the `stretching` argument and the
spacing near the surface is controlled by `refinement`.
"""
function grid_stretching(Lz::Number, Nz::Number, refinement::Number, stretching::Number)

    # Normalise height ranging from 0 to 1
    h(k) = (k - 1) / Nz
    # Linear near-surface generator
    ζ₀(k) = 1 + (h(k) - 1) / refinement
    # Bottom-intensified stretching function
    Σ(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))
    # Generating function
    z_faces(k) = Lz * (ζ₀(k) * Σ(k) - 1)

    return z_faces

end

"""
    function set_two_layer_initial_conditions(model::Oceananigans.AbstractModel,
                                              S::NamedTuple, Θ::NamedTuple;
                                              interface_location = 0.5,
                                              interface_thickness = 100,
                                              salinity_perturbation_width = 100)
Set initial conditions for temperature and salinity in a two layer model. Initial values for
the absolute salinity (`S`) and conservative temperature (`Θ`) in each layer must be passed
to the function as a `NamedTuple` in the form `S = (upper = value, lower = value)` as well
as the model where the initial conditions should be set. **Note** by default the interface
is in the middle of the domain, `interface_location = 0.5`, with width of
`interface_thickness = 100`. A salinity perturbation is added to the upper part of the
domain to trigger turbulent mixing. The width of this perturbation is controlled by the
keyword argument `salinity_pertubration_width`.
"""
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel,
                                           S::NamedTuple, Θ::NamedTuple;
                                           interface_location = 0.5,
                                           interface_thickness = 100,
                                           salinity_perturbation_width = 100)

    ΔS = (S.upper - S.lower) / 2
    ΔΘ = (Θ.upper - Θ.lower) / 2

    initial_S_profile(x, y, z) = ΔS * tanh(interface_thickness * (z + interface_location)) +
                                 (S.lower + ΔS) + perturb_salintiy(z, interface_location,
                                                                salinity_perturbation_width)
    initial_Θ_profile(x, y, z) = ΔΘ * tanh(interface_thickness * (z + interface_location)) +
                                 (Θ.lower + ΔΘ)

    set!(model, S = initial_S_profile, T = initial_Θ_profile)

    return nothing

end

"""
    function perturb_salintiy(z, interface_location)
Where and what value to add to perturb the salinity initial condition.
"""
function perturb_salintiy(z, interface_location, salinity_perturbation_width)
    if z > -interface_location
        exp(-((z + interface_location) - (interface_location / 2))^2 /
              2*(salinity_perturbation_width)^2) / sqrt(2*π*salinity_perturbation_width^2)
    else
        0
    end
end

"""
    function non_dimensional_numbers(model::Oceananigans.AbstractModel,
                                     S::NamedTuple, Θ::NamedTuple;)
Compute and important non-dimensional numbers related to the DNS experiments. The
non-dimensional numbers are:
- Prandtl number: ``Pr = ν / κₜ``
- Schmidt number: ``Sc = ν / κₛ``
- Lewis number:   ``Le = κₜ / κₛ``
- Raleigh number (density): ``Ra_{d} = Ra_{t} / Ra_{s} = (αΔΘ / βΔS) * (1 / Le)``
"""
function non_dimensional_numbers(model::Oceananigans.AbstractModel,
                                 S::NamedTuple, Θ::NamedTuple)

    ν = model.closure.ν
    κₛ, κₜ = model.closure.κ
    Pr = ν / κₜ
    Sc = ν / κₛ
    Le = κₜ / κₛ
    ΔS, ΔΘ = S.upper - S.lower, Θ.upper - Θ.lower
    α, β = gsw_alpha(S.lower, Θ.lower, 0), gsw_beta(S.lower, Θ.lower, 0)
    Ra = ((α * ΔΘ )/ (β * ΔS)) * (1 / Le)

    return Dict("Pr" => Pr, "Sc" => Sc, "Le" => Le, "Ra_ρ" => Ra)

end

"""
    function simulation_progress(sim)
Useful progress messaging for simulation runs
"""
simulation_progress(sim) = @printf("i: % 6d, sim time: % 1.3f, wall time: % 10s, Δt: % 1.4f, advective CFL: %.2e, diffusive CFL: %.2e\n",
                                    iteration(sim), time(sim), prettytime(sim.run_wall_time),
                                    sim.Δt, AdvectiveCFL(sim.Δt)(sim.model),
                                    DiffusiveCFL(sim.Δt)(sim.model))
