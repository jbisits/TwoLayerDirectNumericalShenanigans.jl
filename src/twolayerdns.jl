"""
    struct TwoLayerDNS
Container for all the elements of a `TwoLayerDNS`.
"""
struct TwoLayerDNS{NHM <: NonhydrostaticModel,
                   APF <: AbstractProfileFunction,
                  TLIC <: TwoLayerInitialConditions,
                   ATP <: Union{AbstractTracerPerturbation, Nothing},
                    AN <: Union{AbstractNoise, Nothing}} <: AbstractTwoLayerDNS
    "An [Oceananigans.jl `NonhydrostaticModel`](https://clima.github.io/OceananigansDocumentation/dev/appendix/library/#Oceananigans.Models.NonhydrostaticModels.NonhydrostaticModel-Tuple{})"
                  model :: NHM
    "Continuous profile function"
       profile_function :: APF
    "The two layer initial conditions"
     initial_conditions :: TLIC
    "Perturbation to a tracer field"
    tracer_perturbation :: ATP
    "Initial noise to create instability"
          initial_noise :: AN
end
function Base.show(io::IO, tldns::TwoLayerDNS)
    println(io, "TwoLayerDirectNumericalSimulation")
    println(io, "┣━━━━━━━━━━━━━━━━ model: $(summary(tldns.model))")
    println(io, "┣━━━━━ profile_function: $(typeof(tldns.profile_function))")
    println(io, "┣━━━ initial_conditions: $(typeof(tldns.initial_conditions))")
    println(io, "┣━━ tracer_perturbation: $(typeof(tldns.tracer_perturbation))")
    print(io,   "┗━━━━━━━━ initial_noise: $(typeof(tldns.initial_noise))")
end
TwoLayerDNS(model, profile_function, initial_condition; tracer_perturbation = nothing,
            initial_noise = nothing) =
    TwoLayerDNS(model, profile_function, initial_condition, tracer_perturbation, initial_noise)
Base.iterate(tldns::TwoLayerDNS, state = 1) =
    state > length(fieldnames(TwoLayerDNS)) ? nothing :
                                            (getfield(tldns, state), state + 1)
"""
    function DNSModel(architecture, domain_extent::NamedTuple, resolution::NamedTuple,
                      diffusivities::NamedTuple)
Setup a Direct Numerical Simulation [`Model`](https://clima.github.io/OceananigansDocumentation/dev/appendix/library/#Oceananigans.Models.NonhydrostaticModels.NonhydrostaticModel-Tuple{})
on `architecture` (`CPU()` or `GPU()`) over the
`domain_extent` with `resolution` and scalar diffusivities for momentum
(kinematic viscosity `ν`) and the temperature and salinity tracers (`κₜ` and `κₛ`).
The salinity and temperature tracerse will then be evolved with the equation of state `eos`
chosen.

## Function arguments:

- `architecture`, `CPU()` or `GPU()`;
- `domain_extent` as a `NamedTuple` in the format `(Lx = , Ly = , Lz = )`;
- `resolution` as a `NamedTuple` in the format `(Nx = , Ny = , Nz = )`;
- `diffusivities` as a `NamedTuple` in the format `(ν = , κ = )`, **note** to set different
diffusivities for temperature and salinity `κ` must also be a `NamedTuple` in the format
`κ = (S = , T = )`;
- `eos` a `BoussinesqEquationOfState`, default is [`TEOS10EquationOfState`](https://clima.github.io/SeawaterPolynomials.jl/dev/API/#SeawaterPolynomials.TEOS10.TEOS10EquationOfState)
but any of the [`RoquetEquationOfState`s](https://clima.github.io/SeawaterPolynomials.jl/dev/API/#SeawaterPolynomials.SecondOrderSeawaterPolynomials.RoquetEquationOfState)
may be used.


## Keyword arguments:

- `zgrid_stretching` stretch the grid in the `z` dimension at the bottom of domain at
the rate `stretching`, `false` by default;
- `refinement = 1.2` spacing near the surface in the `z` dimension;
- `stretching = 100` rate of stretching at the bottom of grid in the `z` dimension.
"""
function DNSModel(architecture, domain_extent::NamedTuple, resolution::NamedTuple,
                  diffusivities::NamedTuple, eos::BoussinesqEquationOfState=TEOS10EquationOfState();
                  zgrid_stretching = false,
                  refinement = 1.05,
                  stretching = 40)

    Lx, Ly, Lz = domain_extent.Lx, domain_extent.Ly, domain_extent.Lz
    Nx, Ny, Nz = resolution.Nx, resolution.Ny, resolution.Nz
    zgrid = zgrid_stretching ? grid_stretching(-Lz, Nz, refinement, stretching) : (Lz, 0)

    grid = RectilinearGrid(architecture,
                           topology = (Periodic, Periodic, Bounded),
                           size = (Nx, Ny, Nz),
                           x = (-Lx/2, Lx/2),
                           y = (-Ly/2, Ly/2),
                           z = zgrid)

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
    function TLDNS_simulation_setup(tldns::TwoLayerDNS, Δt::Number, stop_time::Number,
                                    save_schedule::Number;  cfl = 0.75, diffusive_cfl = 0.75,
                                    max_change = 1.2, max_Δt = 1e-1)
Setup a `Simulation` of `tldns` from `initial_conditions` that are of type `TwoLayerInitialConditions`.
Important non-dimensional numbers that are part of this experiment are computed and saved
to the simulation output file.

## Function arguments:

- `model` which is assumed to be a Direct Numerical Simulation setup using `DNS`;
- `Δt` timestep. A timestep wizard is also setup so the size of the timestep may change over
the course of a simulation;
- `stop_time` length of simulation time (in seconds) to run the model for;
- `output_dir` name of the file to save the data to,
- `save_schedule` number (representing time in seconds) at which to save model output, e.g.,
`save_schedule = 1` saves output every second;
- `save_file` a `Symbol` (either `:netcdf` or `:jld2`) choosing whether to save data in
`NetCDF` format (`.nc`) ot `JLD2` format (`.jld2`);
- `save_custom_output!` a function to save custom output during a `Simulation`. **Note:**
the function must be in the form `my_outputs!(simulation, tldns, save_info...; kwargs)`
though can be named anything. See [`save_computed_output!`](@ref) below for an example of how to
write a function to save custom output.

## Keyword arguments:

- `output_path` path where to save output, defualt is `SIMULATION_PATH`
- `checkpointer_time_interval`, pass a `Number` to setup a `Checkpointer` for saving
and restarting a model state at `TimeInterval(Number)`. Defaults to `nothing` in which case
there is no `Checkpointer`.
- `cfl` maximum cfl value used to determine the adaptive timestep size;
- `diffusive_cfl` maximum diffusive cfl value used to determine the adaptive timestep size;
- `max_change` maximum change in the timestep size;
- `max_Δt` the maximum timestep;
- `save_velocities` defaults to `false`, if `true` model velocities will be saved to output;
- `overwrite_saved_output` whether the output overwrites a file of the same name, default is
`true`.
"""
function TLDNS_simulation_setup(tldns::TwoLayerDNS, Δt::Number,
                                stop_time::Number, save_schedule::Number,
                                save_custom_output!::Function=no_custom_output!,
                                save_velocities!::Function=no_velocities!;
                                save_file = :netcdf,
                                output_path = SIMULATION_PATH,
                                checkpointer_time_interval = nothing,
                                cfl = 0.75,
                                diffusive_cfl = 0.75,
                                max_change = 1.2,
                                max_Δt = 1e-1,
                                overwrite_saved_output = true)

    model = tldns.model
    simulation = Simulation(model; Δt, stop_time)

    # time step adjustments
    wizard = TimeStepWizard(; cfl, diffusive_cfl, max_change, max_Δt)
    simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

    # progress reporting
    simulation.callbacks[:progress] = Callback(simulation_progress, IterationInterval(100))

    output_dir = output_directory(tldns, stop_time, output_path)
    save_info = (save_schedule, save_file, output_dir, overwrite_saved_output)

    # model tracers
    save_tracers!(simulation, model, save_info)

    # model velocities
    save_velocities!(simulation, model, save_info)

    # Custom saved output
    save_custom_output!(simulation, tldns, save_info)

    # Checkpointer setup
    checkpointer_setup!(simulation, model, output_dir, checkpointer_time_interval)

    return simulation

end
"""
    function output_directory(tldns::TwoLayerDNS, stop_time::Number, output_path)
Create an `output_directory` for saved output based on the `profile_function`,`initial_conditions`,
`tracer_perturbation` and length of the simulation.
"""
function output_directory(tldns::TwoLayerDNS, stop_time::Number, output_path)

    pf_string = lowercase(string(typeof(tldns.profile_function))[1:findfirst('{', string(typeof(tldns.profile_function))) - 1])
    ic_type = typeof(tldns.initial_conditions)
    ic_string = ic_type <: StableTwoLayerInitialConditions ? "stable" :
                            ic_type <: CabbelingTwoLayerInitialConditions ?
                                "cabbeling" : ic_type <: UnstableTwoLayerInitialConditions ?
                                              "unstable" : ic_type <: IsohalineTwoLayerInitialConditions ?
                                                            "isohaline" : "isothermal"

    tp_string = lowercase(string(typeof(tldns.tracer_perturbation)))
    tp_find = isnothing(findfirst('{', tp_string)) ? length(tp_string) :
                                                     findfirst('{', tp_string) - 1
    stop_time_min = stop_time / 60 ≥ 1 ? string(round(Int, stop_time / 60)) :
                                         string(round(stop_time / 60; digits = 2))
    expt_dir = ic_string *"_"* pf_string *"_"* tp_string[1:tp_find] *"_"* stop_time_min * "min"

    output_dir = mkpath(joinpath(output_path, expt_dir))
    # make a simulation directory if one is not present
    if !isdir(output_dir)
        mkdir(output_dir)
    end

    return output_dir

end
"""
    function save_tracers!(simulation, model, save_schedule, save_file, output_dir, overwrite_saved_output)
Save `model.tracers` during a `Simulation` using an `OutputWriter`.
"""
function save_tracers!(simulation, model, save_schedule, save_file, output_dir,
                       overwrite_saved_output)

    S, T = model.tracers.S, model.tracers.T
    tracers = Dict("S" => S, "T" => T)

    simulation.output_writers[:tracers] =
        save_file == :netcdf ? NetCDFOutputWriter(model, tracers;
                                                  filename = "tracers",
                                                  dir = output_dir,
                                                  overwrite_existing = overwrite_saved_output,
                                                  schedule = TimeInterval(save_schedule)
                                                  ) :
                                JLD2OutputWriter(model, tracers;
                                                 filename = "tracers",
                                                 dir = output_dir,
                                                 schedule = TimeInterval(save_schedule),
                                                 overwrite_existing = overwrite_saved_output)

    return nothing

end
save_tracers!(simulation, model, save_info::Tuple) =
    save_tracers!(simulation, model, save_info...)
"""
    function save_all_velocities!(simulation, model, save_schedule, save_file, output_dir,
                              overwrite_saved_output)
Save `model.velocities` during a `Simulation` using an `OutputWriter`.
"""
function save_all_velocities!(simulation, model, save_schedule, save_file, output_dir,
                          overwrite_saved_output)

    u, v, w = model.velocities
    velocities = Dict("u" => u, "v" => v, "w" => w)

    simulation.output_writers[:velocities] =
        save_file == :netcdf ? NetCDFOutputWriter(model, velocities;
                                                filename = "velocities",
                                                dir = output_dir,
                                                overwrite_existing = overwrite_saved_output,
                                                schedule = TimeInterval(save_schedule)
                                                ) :
                                JLD2OutputWriter(model, velocities;
                                                filename = "velocities",
                                                dir = output_dir,
                                                schedule = TimeInterval(save_schedule),
                                                overwrite_existing = overwrite_saved_output)

    return nothing

end
save_all_velocities!(simulation, model, save_info::Tuple) =
    save_all_velocities!(simulation, model, save_info...)
"""
    function save_vertical_velocities!(simulation, model, save_schedule, save_file, output_dir,
                                    overwrite_saved_output)
Only save vertical velocity.
"""
function save_vertical_velocities!(simulation, model, save_schedule, save_file, output_dir,
                                    overwrite_saved_output)

    w = model.velocities.w
    velocities = Dict("w" => w)

    simulation.output_writers[:velocities] =
    save_file == :netcdf ? NetCDFOutputWriter(model, velocities;
                                filename = "velocities",
                                dir = output_dir,
                                overwrite_existing = overwrite_saved_output,
                                schedule = TimeInterval(save_schedule)
                                ) :
                JLD2OutputWriter(model, velocities;
                                filename = "velocities",
                                dir = output_dir,
                                schedule = TimeInterval(save_schedule),
                                overwrite_existing = overwrite_saved_output)

    return nothing

end
save_vertical_velocities!(simulation, model, save_info::Tuple) =
    save_vertical_velocities!(simulation, model, save_info...)
"Default"
no_velocities!(simulation, model, save_info...) = nothing
"""
    function save_computed_output!(simulation, tldns, save_schedule, save_file, output_dir,
                                   overwrite_saved_output, reference_gp_height)
Save selection of computed output during a `Simulation` using an `OutputWriter`.
"""
function save_computed_output!(simulation, tldns, save_schedule, save_file, output_dir,
                               overwrite_saved_output, reference_gp_height)

    model = tldns.model
    σ = seawater_density(model, geopotential_height = reference_gp_height)
    ϵ = KineticEnergyDissipationRate(model)
    ∫ϵ = Integral(ϵ)
    ϵ_maximum = Reduction(maximum!, ϵ, dims = (1, 2, 3))
    Eₖ = KineticEnergy(model)
    ∫Eₖ = Integral(Eₖ)
    Eₚ = PotentialEnergy(model)
    ∫Eₚ = Integral(Eₚ)

    computed_outputs = Dict("σ" => σ, "∫ϵ" => ∫ϵ, "ϵ_maximum" => ϵ_maximum, "∫Eₖ" => ∫Eₖ,
                            "∫Eₚ" => ∫Eₚ)

    oa = Dict(
        "σ" => Dict("longname" => "Seawater potential density calculated using TEOS-10 at $(reference_gp_height)dbar",
                    "units" => "kgm⁻³"),
        "ϵ_maximum" => Dict("longname" => "Maximum (in space) TKE dissipation"),
        "∫ϵ" => Dict("longname" => "Volume integrated turbulent kintetic energy dissipation"),
        "∫Eₖ" => Dict("longname" => "Volume integrated turbulent kinetic energy"),
        "∫Eₚ" => Dict("longname" => "Volume integrated potential energy (g∫ᵥρzdV)")
        )
    simulation.output_writers[:computed_output] =
        save_file == :netcdf ? NetCDFOutputWriter(model, computed_outputs;
                                                filename = "computed_output",
                                                dir = output_dir,
                                                overwrite_existing = overwrite_saved_output,
                                                schedule = TimeInterval(save_schedule),
                                                output_attributes = oa
                                                ) :
                                JLD2OutputWriter(model, computed_outputs;
                                                filename = "computed_output",
                                                dir = output_dir,
                                                schedule = TimeInterval(save_schedule),
                                                overwrite_existing = overwrite_saved_output)

    non_dimensional_numbers!(simulation, tldns)
    predicted_maximum_density!(simulation, tldns)

    return nothing

end
save_computed_output!(simulation, tldns, save_info::Tuple; reference_gp_height = 0) =
    save_computed_output!(simulation, tldns, save_info..., reference_gp_height)
"Default function for `save_custom_output!` in `TLDNS_simulation_setup`."
no_custom_output!(simulation, model, save_info...) = nothing
"""
    function checkpointer_setup!(simulation, model, output_path, checkpointer_time_interval)
Setup a `Checkpointer` at `checkpointer_time_interval` for a `simulation`
"""
function checkpointer_setup!(simulation, model, output_dir,
                             checkpointer_time_interval::Number)

    checkpoint_dir = joinpath(output_dir, "model_checkpoints/")
    isdir(checkpoint_dir) ? nothing : mkdir(checkpoint_dir)
    schedule = TimeInterval(checkpointer_time_interval)
    cleanup = true
    checkpointer = Checkpointer(model; schedule, dir = checkpoint_dir, cleanup)
    simulation.output_writers[:checkpointer] = checkpointer

    return nothing

end
checkpointer_setup!(simulation, model, output_dir, checkpointer_time_interval::Nothing) = nothing
"""
    function simulation_progress(sim)
Useful progress messaging for simulation runs. This function is from an
[Oceananigans.jl example](https://clima.github.io/OceananigansDocumentation/dev/generated/horizontal_convection/#A-progress-messenger).
"""
simulation_progress(sim) = @printf("i: % 6d, sim time: % 1.3f, wall time: % 10s, Δt: % 1.4f, advective CFL: %.2e, diffusive CFL: %.2e\n",
                                    iteration(sim), time(sim), prettytime(sim.run_wall_time),
                                    sim.Δt, AdvectiveCFL(sim.Δt)(sim.model),
                                    DiffusiveCFL(sim.Δt)(sim.model))
