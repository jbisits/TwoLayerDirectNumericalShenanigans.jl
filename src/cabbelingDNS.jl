"""
    function quasiDNS_cabbeling(resolution::Tuple, diffusivities::NamedTuple)
Setup a quasi (i.e. not true DNS resolution and parameterised `ScalarDiffusivity`) Direct
Numerical Simulation in a 0.5m × 0.5m × 1m box. The argument `resoltuion` is a `NamedTuple`
and needs to be passed as (Nx = , Ny = , Nz = ) and `diffusivities` is a also a `NamedTuple`
of form `(ν = , κ = )`. The equation of state used to evolve the buoyancy is `:Cabbeling`,
from [SeaWaterpolynomials.jl](https://clima.github.io/SeawaterPolynomials.jl/dev/). This
model is designed for experimentation/proof-of-concepts before running a DNS at higher
resolution.
"""
function quasiDNS_cabbeling(resolution::NamedTuple, diffusivities::NamedTuple)

    Lx, Ly, Lz = 0.5, 0.5, 1
    grid = RectilinearGrid(topology = (Periodic, Periodic, Bounded),
                           size = (resolution.Nx, resolution.Ny, resolution.Nz),
                           x = (-Lx/2, Lx/2),
                           y = (-Ly/2, Ly/2),
                           z = (-Lz, 0))

    eos = SeawaterPolynomials.RoquetEquationOfState(:Cabbeling)
    buoyancy = SeawaterBuoyancy(equation_of_state = eos)
    tracers = (:T, :S)

    closure = ScalarDiffusivity(; ν = diffusivities.ν, κ = diffusivities.κ)

    timestepper = :RungeKutta3

    advection = UpwindBiasedFifthOrder()

    return NonhydrostaticModel(; grid, buoyancy, tracers, closure, timestepper, advection)

end
