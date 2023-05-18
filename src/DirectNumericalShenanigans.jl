module DirectNumericalShenanigans

using Oceananigans, SeawaterPolynomials, CairoMakie, JLD2, Printf, Reexport

@reexport using Oceananigans, SeawaterPolynomials, CairoMakie, JLD2, Printf

export SIMULATION_PATH,
       set_two_layer_initial_conditions!

const SIMULATION_PATH = joinpath(@__DIR__, "../data/simulations")

include("twolayermodelsetup.jl")

end
