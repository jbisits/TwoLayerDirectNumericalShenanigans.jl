module TLDNSMakieExt

using CairoMakie, TwoLayerDirectNumericalShenanigans, Printf, NCDatasets,
      GibbsSeaWater, StatsBase
import TwoLayerDirectNumericalShenanigans.erf_tracer_solution
using Oceananigans.Fields

include("plotrecipes.jl")

end
