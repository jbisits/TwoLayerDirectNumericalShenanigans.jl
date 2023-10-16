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
             refinement = 1.05,
             stretching = 40)

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
                                                            TEOS10EquationOfState(;
                                                                reference_density)
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
    function DNS_simulation_setup(dns::TwoLayerDNS, Δt::Number, stop_time::Number,
                                  save_schedule::Number;  cfl = 0.75, diffusive_cfl = 0.75,
                                  max_change = 1.2, max_Δt = 1e-1)
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
`save_schedule = 1` saves output every second;
- `output_writer` a `Symbol` (either `:netcdf` or `:jld2`) choosing whether to save data in
`NetCDF` format (`.nc`) ot `JLD2` format ('.jld2).

## Keyword arguments:

- `cfl` maximum cfl value used to determine the adaptive timestep size;
- `diffusive_cfl` maximum diffusive cfl value used to determine the adaptive timestep size;
- `max_change` maximum change in the timestep size;
- `max_Δt` the maximum timestep;
- `density_reference_gp_height` for the seawater density calculation;
- `save_velocities` defaults to `false`, if `true` model velocities will be saved to output.
"""
function DNS_simulation_setup(dns::TwoLayerDNS, Δt::Number,
                              stop_time::Number, save_schedule::Number,
                              output_writer::Symbol=:netcdf;
                              cfl = 0.75,
                              diffusive_cfl = 0.75,
                              max_change = 1.2,
                              max_Δt = 1e-1,
                              density_reference_gp_height = 0,
                              save_velocities = false)

    model = dns.model
    simulation = Simulation(model; Δt, stop_time)

    # time step adjustments
    wizard = TimeStepWizard(; cfl, diffusive_cfl, max_change, max_Δt)
    simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

    # model tracers
    S, T = model.tracers.S, model.tracers.T

    # Custom saved output
    # Potential density
    parameters = (Zᵣ = density_reference_gp_height,)
    σ = PotentialDensity(model, parameters)

    # Inferred vertical diffusivity
    b_flux = vertical_buoyancy_flux(model)
    ∫ₐb_flux = Integral(b_flux, dims = (1, 2))
    b_grad = ∂b∂z(model)
    ∫ₐb_grad = Integral(b_grad, dims = (1, 2))

    ϵ = KineticEnergyDissipationRate(model)
    # Volume integrated TKE dissipation
    ∫ϵ = Integral(ϵ)
    # Minimum in space Kolmogorov length scale
    η_space(model) = (model.closure.ν^3 / maximum(ϵ))^(1/4)

    # Volume integrated TKE dissipation
    ∫ϵ = Integral(ϵ)

    # Dimensions and attributes for custom saved output
    dims = Dict("η_space" => ())
    oa = Dict(
        "σ" => Dict("longname" => "Seawater potential density calculated using TEOS-10 at $(density_reference_gp_height)dbar",
                    "units" => "kgm⁻³"),
        "η_space" => Dict("longname" => "Minimum (in space) Kolmogorov length"),
        "∫ₐb_flux" => Dict("longname" => "Horizontally integrated vertical buoyacny flux (w'b')"),
        "∫ₐb_grad" => Dict("longname" => "Horizontally integrated vertical buoyancy gradient (∂b/∂z)"),
        "∫ϵ" => Dict("longname" => "Volume integrated turbulent kintetic energy dissipation")
        )

    # outputs to be saved during the simulation
    outputs = Dict("S" => S, "T" => T, "η_space" => η_space, "σ" => σ, "∫ϵ" => ∫ϵ,
                  "∫ₐb_flux" => ∫ₐb_flux, "∫ₐb_grad" => ∫ₐb_grad)
    if save_velocities
        u, v, w = model.velocities
        velocities = Dict("u" => u, "v" => v, "w" => w)
        merge!(outputs, velocities)
    end

    filename = form_filename(dns, stop_time, output_writer)
    simulation.output_writers[:outputs] = output_writer == :netcdf ?
                                            NetCDFOutputWriter(model, outputs,
                                                            filename = filename,
                                                            schedule = TimeInterval(save_schedule),
                                                            overwrite_existing = true,
                                                            dimensions = dims,
                                                            output_attributes = oa
                                                            ) :
                                            JLD2OutputWriter(model, outputs,
                                                            filename = filename,
                                                            schedule = TimeInterval(save_schedule),
                                                            overwrite_existing = true)

    non_dimensional_numbers!(simulation, dns)
    predicted_maximum_density!(simulation, dns)

    # progress reporting
    simulation.callbacks[:progress] = Callback(simulation_progress, IterationInterval(100))

    return simulation

end
"""
    function form_filename(dns::TwoLayerDNS, stop_time::Number, output_writer::Symbol)
Create a filename for saved output based on the `profile_function`,`initial_conditions`,
`tracer_perturbation` and length of the simulation.
"""
function form_filename(dns::TwoLayerDNS, stop_time::Number, output_writer::Symbol)

    pf_string = lowercase(string(typeof(dns.profile_function))[1:findfirst('{', string(typeof(dns.profile_function))) - 1])
    ic_type = typeof(dns.initial_conditions)
    ic_string = ic_type <: StableTwoLayerInitialConditions ? "stable" :
                            ic_type <: CabbelingTwoLayerInitialConditions ?
                                "cabbeling" : ic_type <: UnstableTwoLayerInitialConditions ?
                                              "unstable" : ic_type <: IsohalineTwoLayerInitialConditions ?
                                                            "isohaline" : "isothermal"

    tp_string = lowercase(string(typeof(dns.tracer_perturbation)))
    tp_find = isnothing(findfirst('{', tp_string)) ? length(tp_string) :
                                                     findfirst('{', tp_string) - 1
    stop_time_min = stop_time / 60 ≥ 1 ? string(round(Int, stop_time / 60)) :
                                         string(round(stop_time / 60; digits = 2))
    filetype = output_writer == :netcdf ? ".nc" : ".jld2"
    savefile = ic_string *"_"* pf_string *"_"* tp_string[1:tp_find] *"_"* stop_time_min * "min" * filetype

    # make a simulation directory if one is not present
    if !isdir(SIMULATION_PATH)
        mkdir(SIMULATION_PATH)
    end
    filename = joinpath(SIMULATION_PATH, savefile)

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
