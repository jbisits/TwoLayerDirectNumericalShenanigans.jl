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

const SIMULATION_PATH = joinpath(pwd(), "data/simulations")
const DNCS = DirectNumericalCabbelingShenanigans # alias

include("DNS.jl")
include("twolayersetup.jl")
include("output_utils.jl")

end
