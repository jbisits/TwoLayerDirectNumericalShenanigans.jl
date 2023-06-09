module DirectNumericalCabbelingShenanigans

using Oceananigans, SeawaterPolynomials, CairoMakie, JLD2, Printf, Reexport, GibbsSeaWater

@reexport using Oceananigans, SeawaterPolynomials, CairoMakie, JLD2, Printf, GibbsSeaWater

export SIMULATION_PATH,
       DNS_cabbeling,
       set_two_layer_initial_conditions!,
       non_dimensional_numbers,
       simulation_progress,
       animate_2D_field

const SIMULATION_PATH = joinpath(@__DIR__, "../data/simulations")

if !isdir(SIMULATION_PATH)
    mkdir(SIMULATION_PATH)
end

include("cabbelingDNS.jl")
include("utils.jl")

end
