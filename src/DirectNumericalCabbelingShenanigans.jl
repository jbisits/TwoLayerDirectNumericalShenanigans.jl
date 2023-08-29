module DirectNumericalCabbelingShenanigans

using Oceananigans, Printf, Reexport
using SeawaterPolynomials: TEOS10EquationOfState
using Oceananigans: AbstractModel

@reexport using Oceananigans, Reexport

export
    SIMULATION_PATH,
    DNCS,
    DNS,
    DNS_simulation_setup

export
    animate_2D_field,
    visualise_initial_conditions,
    visualise_initial_density,
    visualise_snapshot

const SIMULATION_PATH = joinpath(pwd(), "data/simulations")
const DNCS = DirectNumericalCabbelingShenanigans # alias

include("DNS.jl")
include("twolayersetup.jl")

end
