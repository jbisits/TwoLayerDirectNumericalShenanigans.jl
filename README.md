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
julia> using TwoLayerDirectNumericalShenanigans
```

The package then allows the setup of Direct Numerical Simulations with two layers ideal for exploring double diffusion or other small scale processes.
Documentation will come soon!
