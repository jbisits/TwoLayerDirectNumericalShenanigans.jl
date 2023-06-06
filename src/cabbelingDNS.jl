"""
    function DNS_cabbeling(resolution::Tuple, diffusivities::NamedTuple)
Setup a Direct Numerical Simulation in a 0.05m × 0.05m × 1m box.
The argument `resoltuion` is a `NamedTuple` and needs to be passed as (Nx = , Ny = , Nz = )
and `diffusivities` is a also a `NamedTuple` of form `(ν = , κ = )` for the momentum and
tracer diffusivities respectively (salt and heat can have different diffusivities).
The equation of state used to evolve the buoyancy is the [polynomial approximation to the
full non-linear TEOS10 equation of state appropriate for Boussinesq models](https://clima.github.io/SeawaterPolynomials.jl/dev/API/#SeawaterPolynomials.TEOS10).
The keyword argument `reference_density` can be passed to evolve the density
about. If nothing is passed default reference density of 1020kgm⁻³ is used.
"""
function DNS_cabbeling(resolution::NamedTuple, diffusivities::NamedTuple;
                       reference_density = nothing)

    Lx, Ly, Lz = 0.1, 0.1, 1
    Nx, Ny, Nz = resolution.Nx, resolution.Ny, resolution.Nz

    # Grid spacing
    refinement = 1.2 # controls spacing near surface (higher means finer spaced)
    stretching = 100  # controls rate of stretching at bottom

    # Normalized height ranging from 0 to 1
    h(k) = (k - 1) / resolution.Nz
    # Linear near-surface generator
    ζ₀(k) = 1 + (h(k) - 1) / refinement
    # Bottom-intensified stretching function
    Σ(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))
    # Generating function
    z_faces(k) = 1 * (ζ₀(k) * Σ(k) - 1)

    grid = RectilinearGrid(topology = (Periodic, Periodic, Bounded),
                           size = (Nx, Ny, Nz),
                           x = (-Lx/2, Lx/2),
                           y = (-Ly/2, Ly/2),
                           z = z_faces)

    eos = isnothing(reference_density) ? SeawaterPolynomials.TEOS10EquationOfState() :
                            SeawaterPolynomials.TEOS10EquationOfState(; reference_density)
    buoyancy = SeawaterBuoyancy(equation_of_state = eos)
    tracers = (:T, :S)

    closure = ScalarDiffusivity(; ν = diffusivities.ν, κ = diffusivities.κ)

    timestepper = :RungeKutta3

    advection = CenteredSecondOrder()

    return NonhydrostaticModel(; grid, buoyancy, tracers, closure, timestepper, advection)

end
