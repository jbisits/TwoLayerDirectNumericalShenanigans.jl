"""
    function animate_2D_field(field_timeseries::FieldTimeSeries, field_name::AbstractString,
                              field_dimensions::NamedTuple{Symbol, Symbol}; colormap = :thermal)
Animate the time series of `field_name` a `FieldTimeSeries` over `field_dimensions`.
"""
function animate_2D_field end
"""
    function visualise_initial_conditions(dns::TwoLayerDNS, xslice::Integer, yslice::Integer)
Plot the initial state of the `tracers` in a `model`. This function assumes there are two
tracers (salinity and temperature) and plots the x-z, y-z and field-z initial fields at
`xslice` and `yslice`.
"""
function visualise_initial_conditions end
"""
    function visualise_initial_stepchange(dns::TwoLayerDNS, interface_location::Number)
Plot an initial step change of the `tracers` in a `model`. This function assumes there are two
tracers (salinity and temperature) and plots the x-z, y-z and field-z initial fields.
"""
function visualise_initial_stepchange end
"""
    initial_tracer_heaviside(z, C::Number, ΔC::Number, interface_location)
Modified Heaviside function for initial condition of a tracer with depth. The interface_location of the
Heaviside function is calculated from the `extrema` of the depth array `z`.

## Function arguments:

- `z` for the Oceananigans model grid to evaulate the function at;
- `C` tracer value in deeper part of the step;
- `ΔC` difference in tracer between the steps.

## Keyword arguments:
- `interface_location` where the step takes place e.g. `z - interface_location < 0 ? 0 : 1`.
Default behaviour puts the `interface_location` in the centre of the depth range given by `z`.
"""
function initial_tracer_heaviside end
"""
    function visualise_initial_density(dns::TwoLayerDNS, xslice::Integer,  yslice::Integer,
                                       pressure::Union{Number, Vector{Number}})
Compute and plot the initial density at `pressure` (either reference pressure or in-situ
pressure). The arguments `xslice` and `yslice` are used to choose where in the domain the
figures are from.
"""
function visualise_initial_density end
"""
function visualise_snapshot(field_timeseries::FieldTimeSeries, field_name::AbstractString,
                            snapshot::Int64)
Plot a `snapshot` of the `field_timeseries`  with `field_name` at `xslice`, `yslice`.
"""
function visualise_snapshot end
"""
    function animate_tracer_distributions(tracers::AbstractString;
                                          S_binwidth = 0.001, Θ_binwidth = 0.01)
Animate volume distribution of the salinity and temperature in `tracers`.
**Note:** this assumens `tracers` is a `.nc` file.
"""
function animate_tracer_distributions end
"""
    function animate_joint_tracer_distribution(tracers::AbstractString;
                                               S_binwidth = 0.001, Θ_binwidth = 0.01)
Animate the joint distribution of the salinity and temperature in `tracers`.
**Note:** this assumens `tracers` is a `.nc` file.
"""
function animate_joint_tracer_distribution end
"""
    function plot_scalar_diagnostics(computed_output::AbstractString)
Plot scalar diagnostics from `computed_output`. **Note:**
this function assumes that `computed_output` is a `.nc` file and is intended to be extended
based on what scalar diagnostics are computed during the simulation.
"""
function plot_scalar_diagnostics end
"""
    function hovmoller(computed_output::AbstractString, variable::AbstractString)
Produce a Hovmoller plot (with time on xaxis and z on yaxis) of a `variable` in `computed_output`.
"""
function hovmoller end
"""
    function animate_tracers(tracers::AbstractString)
Animate the salinity and temperature `tracers` from saved `.nc` output.
"""
function animate_tracers end
"""
    function animate_density(computed_output::AbstractString, variable::AbstractString;
                                     xslice = 52, yslice = 52)
Animate the density `variable` in `computed_output`. The `xslice` and `yslice` keyword
arguments can be passed to specify where in the domain the slices are taken from.
"""
function animate_density end
