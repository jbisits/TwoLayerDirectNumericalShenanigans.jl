"""
    const DOMAIN_EXTENT
Domain extent on which the two layer simulations are run.
"""
const DOMAIN_EXTENT = (Lx = 0.1, Ly = 0.1, Lz = 1)
"""
    const HIGH_RESOLUTION
Resolution at which to run the DNS sufficient to resolve turbulence on all scales, i.e.
`dx`, `dy`, `dz` < `η_min` where `η_min` is smallest Kolmogorov scale in length and time for
the run of a simulation.
"""
const HIGH_RESOLUTION = (Nx = 124, Ny = 124, Nz = 1400)
"""
    const SO_DIFFUSIVITIES
Diffusivity estimates for the Southern Ocean.
"""
const SO_DIFFUSIVITIES = (ν = 1e-6, κ = (S = 1e-9, T = 1e-7))
"""
    const REFERENCE_DENSITY
Reference density for use in the two layer DNS. Calculated using the salinity `S₀ˡ` and
temperature `T₀ˡ` of the lower layer .
"""
const REFERENCE_DENSITY = gsw_rho(34.7, 0.5, 0)
"""
    const INTERFACE_LOCATION
Location of the interface (in the vertical) between the upper and lower layers.
"""
const INTERFACE_LOCATION = -0.375
"""
    const SIMULATION_PATH
Path to where the simulations are saved by default. If the folder does not exist it will be
created when initialising a `Simulation` with `DNS_simulation_setup`.
"""
const SIMULATION_PATH = joinpath(pwd(), "data/simulations")
"""
    const CHECKPOINT_PATH
Path to where `Checkpoint`s are saved if used in a simulation.
"""
const CHECKPOINT_PATH = joinpath(SIMULATION_PATH, "model_checkpoints/")
"Alias for DirectNumericalCabbelingShenanigans."
const TLDNS = TwoLayerDirectNumericalShenanigans # alias
