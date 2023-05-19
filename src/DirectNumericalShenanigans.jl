module DirectNumericalShenanigans

using Oceananigans, SeawaterPolynomials, CairoMakie, JLD2, Printf, Reexport

@reexport using Oceananigans, SeawaterPolynomials, CairoMakie, JLD2, Printf

export SIMULATION_PATH,
       quasiDNS_cabbeling,
       set_two_layer_initial_conditions!,
       simulation_progress

const SIMULATION_PATH = joinpath(@__DIR__, "../data/simulations")

if !isdir(SIMULATION_PATH)
    mkdir(SIMULATION_PATH)
end

include("cabbelingDNS.jl")
include("twolayermodelsetup.jl")

end
