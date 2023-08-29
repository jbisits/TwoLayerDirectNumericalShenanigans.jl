"""
    struct TwoLayerDNS
Container for all the elements of a `TwoLayerDNS` with a `SalinityPerturbation`.
"""
struct TwoLayerDNS{NHM <: NonhydrostaticModel, CPF <: ContinuousProfileFunction,
                   TLIC <: TwoLayerInitialConditions, SP <: SalinityPerturbation}
    "An [Oceananigans.jl `NonhydrostaticModel`](https://clima.github.io/OceananigansDocumentation/dev/appendix/library/#Oceananigans.Models.NonhydrostaticModels.NonhydrostaticModel-Tuple{})"
    model :: NHM
    "Continuous profile function"
    profile_function :: CPF
    "The two layer initial conditions"
    initial_conditions :: TLIC
    "Perturbation to the salinity"
    salinity_perturbation :: SP
end
function Base.show(io::IO, tldns::TwoLayerDNS)
    println(io, "TwoLayerDirectNumericalSimulation")
    println(io, " ┣━━━━━━━━━━━━━━━━━ model: $(summary(tldns.model))")
    println(io, " ┣━━━━━━ profile_function: $(typeof(tldns.profile_function))")
    println(io, " ┣━━━━ initial_conditions: $(typeof(tldns.initial_conditions))")
    print(io,   " ┗━ salinity_perturbation: $(typeof(tldns.salinity_perturbation))")
end
"""
    function DNS(architecture, domain_extent::NamedTuple, resolution::NamedTuple,
                 diffusivities::NamedTuple)
Setup a Direct Numerical Simulation on `architecture` (`CPU()` or `GPU()`) over the
`domain_extent` with `resolution` and scalar diffusivities for momentum
(kinematic viscosity `ν`) and the temperature and salinity tracers (`κₜ` and `κₛ`).
To evolve the temperature and salinity tracers, the [polynomial approximation to the
full non-linear TEOS10 equation of state appropriate for Boussinesq models](https://clima.github.io/SeawaterPolynomials.jl/dev/API/#SeawaterPolynomials.TEOS10)
is used

## Function arguments:

- architecture, `CPU()` or `GPU()`;
- domain_extent::NamedTuple in the format `(Lx = , Ly = , Lz = )`;
- resolution::NamedTuple in the format `(Nx = , Ny = , Nz = )`;
- diffusivities::NamedTuple in the format `(ν = , κ = )`.

**Note:** to set different diffusivities for temperature and salinity `κ` must also be a
`NamedTuple` in the format `κ = (S = , T = )`.

## Keyword arguments:

- `linear_eos`, set a linear equation of state for the DNS, default is `false`.
- `α`, set thermal expansion coefficient for use with `linear_eos`, default value is the
same as set in [Oceananigans.jl](https://clima.github.io/OceananigansDocumentation/dev/appendix/library/#Oceananigans.BuoyancyModels.LinearEquationOfState).
- `β`, set haline contraction coefficient for use with `linear_eos`, default value is the
same as set in [Oceananigans.jl](https://clima.github.io/OceananigansDocumentation/dev/appendix/library/#Oceananigans.BuoyancyModels.LinearEquationOfState).
- `reference_density = nothing` for use in the full non-linear equation of state,
defaults to 1020kgm⁻³ but any value can be passed in;
- `zgrid_stretching = true` stretch the grid in the `z` dimension at the bottom of domain at
the rate `stretching`, if `false` uniform grid spacing is used;
- `refinement = 1.2` spacing near the surface in the `z` dimension;
- `stretching = 100` rate of stretching at the bottom of grid in the `z` dimension.
"""
function DNS(architecture, domain_extent::NamedTuple, resolution::NamedTuple,
             diffusivities::NamedTuple;
             linear_eos = false,
             α = 1.67e-4,
             β = 7.80e-4,
             reference_density = nothing,
             zgrid_stretching = true,
             refinement = 1.2,
             stretching = 100)

    Lx, Ly, Lz = domain_extent.Lx, domain_extent.Ly, domain_extent.Lz
    Nx, Ny, Nz = resolution.Nx, resolution.Ny, resolution.Nz
    zgrid = zgrid_stretching ? grid_stretching(Lz, Nz, refinement, stretching) : (-Lz, 0)

    grid = RectilinearGrid(architecture,
                           topology = (Periodic, Periodic, Bounded),
                           size = (Nx, Ny, Nz),
                           x = (-Lx/2, Lx/2),
                           y = (-Ly/2, Ly/2),
                           z = zgrid)

    eos = linear_eos == true ? LinearEquationOfState(thermal_expansion = α,
                                                    haline_contraction = β) :
                                                    isnothing(reference_density) ?
                                                            TEOS10EquationOfState() :
                                                            TEOS10EquationOfState(; reference_density)
    buoyancy = SeawaterBuoyancy(equation_of_state = eos)
    tracers = (:S, :T)

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
    function DNS_simulation_setup(model::Oceananigans.AbstractModel, Δt::Number,
                                  stop_time::Number,
                                  initial_conditions::TwoLayerInitialConditions)
Setup a DNS from `initial_conditions` that are of type `TwoLayerInitialConditions`.
Important non-dimensional numnbers that are part of this experiment are computed and saved
to the simulation output file.

## Function arguments:

- `model` which is assumed to be a Direct Numerical Simulation setup using `DNS`;
- `Δt` timestep. A timestep wizard is also setup so the size of the timestep may change over
the course of a simulation;
- `stop_time` length of simulation time (in seconds) to run the model for;
- `savefile` name of the file to save the data to,
- `save_schedule` number (representing time in seconds) at which to save model output, e.g.,
`save_schedule = 1` saves output every second.

## Keyword arguments:

- `cfl` maximum cfl value used to determine the adaptive timestep size;
- `diffusive_cfl` maximum diffusive cfl value used to determine the adaptive timestep size;
- `max_change` maximum change in the timestep size;
- `max_Δt` the maximum timestep.
"""
function DNS_simulation_setup(model::Oceananigans.AbstractModel, Δt::Number,
                              stop_time::Number, save_schedule::Number,
                              initial_conditions::TwoLayerInitialConditions;
                              cfl = 0.75,
                              diffusive_cfl = 0.75,
                              max_change = 1.2,
                              max_Δt = 1e-1)

    simulation = Simulation(model; Δt, stop_time)

    # time step adjustments
    wizard = TimeStepWizard(; cfl, diffusive_cfl, max_change, max_Δt)
    simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

    # save output
    ϵ = KineticEnergyDissipationRate(model)
    outputs = (S = model.tracers.S, T = model.tracers.T, ϵ = ϵ, w = model.velocities.w)
    filename = form_filename(initial_conditions)
    simulation.output_writers[:outputs] = JLD2OutputWriter(model, outputs,
                                                    filename = filename,
                                                    schedule = TimeInterval(save_schedule),
                                                    overwrite_existing = true)
    jldopen(filename, "a+") do file
        file["Non_dimensional_numbers"] = non_dimensional_numbers(model, initial_conditions)
    end

    # progress reporting
    simulation.callbacks[:progress] = Callback(simulation_progress, IterationInterval(100))

    return simulation

end
"""
    function form_filename(initial_conditions::TwoLayerInitialConditions)
Create a directory based on the temperature of the upper layer and a file for the saved
output based on the type of initial condition (i.e. stable, cabbeling or unstable).
"""
function form_filename(initial_conditions::TwoLayerInitialConditions)

    ic_type = typeof(initial_conditions)
    savefile = ic_type <: StableTwoLayerInitialConditions ? "stable" :
                            ic_type <: CabbelingTwoLayerInitialConditions ?
                                "cabbeling" : ic_type <: UnstableTwoLayerInitialConditions ?
                                              "unstable" : ic_type <: IsohalineTwoLayerInitialConditions ?
                                                            "isohaline" : "isothermal"
    # make a simulation directory if one is not present
    if !isdir(SIMULATION_PATH)
        mkdir(SIMULATION_PATH)
    end
    filename = joinpath(SIMULATION_PATH, savefile * ".jld2")

    return filename

end
"""
    function simulation_progress(sim)
Useful progress messaging for simulation runs. This function is from an
[Oceananigans.jl example](https://clima.github.io/OceananigansDocumentation/dev/generated/horizontal_convection/#A-progress-messenger).
"""
simulation_progress(sim) = @printf("i: % 6d, sim time: % 1.3f, wall time: % 10s, Δt: % 1.4f, advective CFL: %.2e, diffusive CFL: %.2e\n",
                                    iteration(sim), time(sim), prettytime(sim.run_wall_time),
                                    sim.Δt, AdvectiveCFL(sim.Δt)(sim.model),
                                    DiffusiveCFL(sim.Δt)(sim.model))
