module DirectNumericalCabbelingShenanigans

using Oceananigans, SeawaterPolynomials, Printf, Reexport
using Oceananigans: AbstractModel

@reexport using Oceananigans, Reexport

export
    SIMULATION_PATH,
    DNCS,
    DNS,
    DNS_simulation_setup

const SIMULATION_PATH = joinpath(@__DIR__, "../data/simulations")
const DNCS = DirectNumericalCabbelingShenanigans # alias

if !isdir(SIMULATION_PATH)
    mkdir(SIMULATION_PATH)
end

include("DNS.jl")
include("twolayersetup.jl")
include("output_utils.jl")

end
