module DirectNumericalShenanigans

using Oceananigans, SeawaterPolynomials, CairoMakie, JLD2, Printf, Reexport

@reexport using Oceananigans, SeawaterPolynomials, CairoMakie, JLD2, Printf

export SIMULATION_PATH

const SIMULATION_PATH = joinpath(@__DIR__, "../data/simulations")

end
