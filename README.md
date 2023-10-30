# TwoLayerDirectNumericalShenanigans.jl

[![Build Status](https://github.com/jbisits/TwoLayerDirectNumericalShenanigans.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jbisits/TwoLayerDirectNumericalShenanigans.jl/actions/workflows/CI.yml?query=branch%3Amain)

A Julia package (unregistered) to explore two layer systems using Direct Numerical Simulation experiments built with [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl).

To add the package (assuming julia is already installed):

```julia
julia> using Pkg
julia> Pkg.add(url="https://github.com/jbisits/TwoLayerDirectNumericalShenanigans.jl")
```

To then use the package

```julia
julia> using TwoLayerDirectNumericalShenanigans
```

The package then allows the setup of Direct Numerical Simulations with two layers ideal for exploring double diffusion or other small scale processes.
Documentation will come soon!
