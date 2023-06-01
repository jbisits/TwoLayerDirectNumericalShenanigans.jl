"""
    function quasiDNS_cabbeling(resolution::Tuple, diffusivities::NamedTuple)
Setup a quasi (i.e. not true DNS resolution and parameterised `ScalarDiffusivity`) Direct
Numerical Simulation in a 0.05m × 0.05m × 1m box. The argument `resoltuion` is a `NamedTuple`
and needs to be passed as (Nx = , Ny = , Nz = ) and `diffusivities` is a also a `NamedTuple`
of form `(ν = , κ = )` for the momentum and tracer diffusivities.
The equation of state used to evolve the buoyancy is the [polynomial approximation to the
full non-linear TEOS10 equation of state appropriate for Boussinesq models](https://clima.github.io/SeawaterPolynomials.jl/dev/API/#SeawaterPolynomials.TEOS10).
This model is designed for experimentation/proof-of-concepts before running a DNS at higher
resolution. The keyword argument `reference_density` can be passed to evolve the density
about. If nothing is passed default reference density of 1020kgm⁻³ is used.
"""
function quasiDNS_cabbeling(resolution::NamedTuple, diffusivities::NamedTuple;
                            reference_density = nothing)

    Lx, Ly, Lz = 0.1, 0.1, 1
    grid = RectilinearGrid(topology = (Periodic, Periodic, Bounded),
                           size = (resolution.Nx, resolution.Ny, resolution.Nz),
                           x = (-Lx/2, Lx/2),
                           y = (-Ly/2, Ly/2),
                           z = (-Lz, 0))

    eos = isnothing(reference_density) ? SeawaterPolynomials.TEOS10EquationOfState() :
                            SeawaterPolynomials.TEOS10EquationOfState(; reference_density)
    buoyancy = SeawaterBuoyancy(equation_of_state = eos)
    tracers = (:T, :S)

    closure = ScalarDiffusivity(; ν = diffusivities.ν, κ = diffusivities.κ)

    timestepper = :RungeKutta3

    advection = WENO()

    return NonhydrostaticModel(; grid, buoyancy, tracers, closure, timestepper, advection)

end
