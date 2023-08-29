"""
    const DOMAIN_EXTENT
Domain extent on which the two layer simulations are run.
"""
const DOMAIN_EXTENT = (Lx = 0.1, Ly = 0.1, Lz = 1)
"""
    const HIGH_RESOLUTION
Resolution (high) at which to run the DNS.
"""
const HIGH_RESOLUTION = (Nx = 50, Ny = 50, Nz = 1400)
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
const REFERENCE_DENSITY = gsw_rho(34.7, -1.5, 0)
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
"Alias for DirectNumericalCabbelingShenanigans."
const DNCS = DirectNumericalCabbelingShenanigans # alias
