# DirectNumericalCabbelingShenanigans.jl

[![Build Status](https://github.com/jbisits/DirectNumericalShenanigans.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jbisits/DirectNumericalShenanigans.jl/actions/workflows/CI.yml?query=branch%3Amain)

A Julia package (unregistered) to explore the cabbeling instability using Direct Numerical Simulation experiments built with [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl).

To add the package (assuming julia is already installed):

```julia
julia> using Pkg
julia> Pkg.add(url="https://github.com/jbisits/DirectNumericalCabbelingShenanigans.jl.git")
```

To then use the package

```julia
julia> using DirectNumericalCabbelingShenanigans
```

The main DNS setup will then be able to be used.
There are two sub modules (`TwoLayerDNS` and `OutputAnalysis`) that are also available but need to be called explicitly if you wish to use them.

## `TwoLayerDNS`

Utilities to help setup two layer DNS experiments.

## `OutputAnalysis`

Some functions for plotting model output.
