"""
    module OutputUtilities
Module to help process and analyse output data from running a Direct Numerical Simulation.
"""
module OutputUtilities

using DirectNumericalCabbelingShenanigans, Printf
using Oceananigans.Fields

@reexport using CairoMakie, JLD2, GibbsSeaWater

export
    animate_2D_field,
    visualise_initial_conditions,
    visualise_initial_density,
    visualise_snapshot,
    compute_density

"""
    function animate_2D_field(field_timeseries::FieldTimeSeries, field_name::AbstractString,
                              field_dimensions::NamedTuple{Symbol, Symbol}; colormap = :thermal)
Animate a time series that is saved in memory.

Function arguments:

- `field_timeseries` to be animated;
- `field_name` the name of the variable in the `field_timeseries`
- `field_dimensions` the dimensions over which to animate, two dimensions need to be passed
in a `Tuple` where each entry is a `Symbol` e.g. `(:x, :z)`.

Keyword arguments:

- `colormap` for the animated `field_timeseries`
"""
function animate_2D_field(field_timeseries::FieldTimeSeries, field_name::AbstractString,
                          field_dimensions::Tuple{Symbol, Symbol}; colormap = :thermal)

    x, y, z = nodes(field_timeseries[1])
    plot_dims = field_dimensions == (:x, :z) ? (x, z) :
                                               field_dimensions == (:y, :z) ? (y, z) : (x, y)
    t = field_timeseries.times
    n = Observable(1)
    field_tₙ = @lift interior(field_timeseries[$n], :, 1, :)
    c_limits = extrema(interior(field_timeseries, :, :, :, 1))
    title = @lift @sprintf("t=%1.2f", t[$n])
    fig, ax, hm = heatmap(plot_dims[1], plot_dims[2], field_tₙ;
                          colormap, colorrange = c_limits)
    ax.xlabel = string(field_dimensions[1])
    ax.ylabel = string(field_dimensions[2])
    Colorbar(fig[1, 2], hm, label = field_name)

    frames = eachindex(t)
    filename = string(field_dimensions[1]) * string(field_dimensions[2]) * "_" *
                field_name
    record(fig, joinpath(@__DIR__, "../data/analysis/", filename * ".mp4"),
          frames, framerate=8) do i
        msg = string("Plotting frame ", i, " of ", frames[end])
        print(msg * " \r")
        n[] = i
    end

end

"""
    function visualise_initial_conditions(model::Oceanangians.AbstractModel)
Plot the initial state of the `tracers` in a `model`. This function assumes there are two
tracers (salinity and temperature) and plots the x-z, y-z and field-z initial fields.
"""
function visualise_initial_conditions(model::Oceananigans.AbstractModel)

    x, y, z = nodes(model.grid, (Center(), Center(), Center()))
    S = model.tracers.S
    T = model.tracers.T
    fig = Figure(size = (600, 1500))
    ax = [Axis(fig[j, i]) for i ∈ 1:2, j ∈ 1:3]

    hm = heatmap!(ax[1], x, z, interior(S, :, 1, :, 1); colormap = :haline)
    ax[1].title = "Initial salinity (x-z)"
    ax[1].xlabel = "x (m)"
    ax[1].ylabel = "z (m)"
    heatmap!(ax[2], x, z, interior(S, 1, :, :, 1); colormap = :haline)
    ax[2].title = "Initial salinity (y-z)"
    ax[2].xlabel = "y (m)"
    ax[2].ylabel = "z (m)"
    Colorbar(fig[1, 3], hm, label = "S (gkg⁻¹)")
    hm = heatmap!(ax[3], y, z, interior(T, :, 1, :, 1); colormap = :thermal)
    ax[3].title = "Initial temperature (x-z)"
    ax[3].xlabel = "x (m)"
    ax[3].ylabel = "z (m)"
    heatmap!(ax[4], y, z, interior(T, 1, :, :, 1); colormap = :thermal)
    ax[4].title = "Initial temperature (y-z)"
    ax[4].xlabel = "y (m)"
    ax[4].ylabel = "z (m)"
    Colorbar(fig[2, 3], hm, label = "Θ (°C)")
    lines!(ax[5], interior(S, 1, 1, :, 1), z)
    ax[5].title = "Initial salinity profile"
    ax[5].xlabel = "S (gkg⁻¹)"
    ax[5].ylabel = "z (m)"
    lines!(ax[6], interior(T, 1, 1, :, 1), z)
    ax[6].title = "Initial temperature profile"
    ax[6].xlabel =  "Θ (°C)"
    ax[6].ylabel = "z (m)"

    return fig

end
"""
    function visualise_initial_density(model::Oceananigans.AbstractModel,
                                       pressure::Union{Number, Vector{Number}})
Compute and plot the initial density against depth at `pressure` (either reference pressure
or in-situ pressure).
"""
function visualise_initial_density(model::Oceananigans.AbstractModel,
                                   pressure::Union{Number, Vector{Number}})

    x = xnodes(model.grid, Center(), Center(), Center())
    z = znodes(model.grid, Center(), Center(), Center())
    S = interior(model.tracers.S, :, 1, :, 1)
    T = interior(model.tracers.T, :, 1, :, 1)
    ρ = gsw_rho.(S, T, pressure)

    fig = Figure(size = (1000, 600))
    ax = [Axis(fig[1, i]) for i ∈ 1:2]

    hm = heatmap!(ax[1], x, z, ρ; colormap = :dense)
    ax[1].title = "Initial density (x-z)"
    ax[1].xlabel = "x (m)"
    ax[1].ylabel = "z (m)"
    Colorbar(fig[2, 1], hm, label = "ρ (kgm⁻³)", vertical = false, flipaxis = false)
    lines!(ax[2], ρ[1, :], z)
    ax[2].title = "Initial density profile"
    ax[2].xlabel = "ρ (kgm⁻³)"
    ax[2].ylabel = "z (m)"

    return fig

end

"""
function visualise_snapshot(field_timeseries::FieldTimeSeries, field_name::AbstractString,
                            snapshot::Int64)
Plot a `snapshot` of the `field_timeseries`  with `field_name`.
"""
function visualise_snapshot(field_timeseries::FieldTimeSeries, field_name::AbstractString,
                            snapshot::Int64; colormap = :thermal)

    x, y, z = nodes(field_timeseries[1])
    t = round(field_timeseries.times[snapshot]; digits = 3)
    fig = Figure(size = (1000, 500))
    ax = [Axis(fig[1, i]) for i ∈ 1:2]

    hm = heatmap!(ax[1], x, z, interior(field_timeseries, :, 1, :, snapshot); colormap)
    ax[1].title = field_name * " at time t = $(t) (x-z)"
    ax[1].xlabel = "x (m)"
    ax[1].ylabel = "z (m)"
    Colorbar(fig[2, 1], hm, vertical = false, label = field_name, flipaxis = false)
    lines!(ax[2], interior(field_timeseries, 1, 1, :, snapshot), z)
    ax[2].title = field_name * " profile at time t = $(t)"
    ax[2].xlabel = field_name
    ax[2].ylabel = "z (m)"

    return fig

end

"""
    function compute_density(S_timeseries::FieldTimeSeries, T_timeseries::FieldTimeSeries;
                             reference_pressure = 0)
Return a density `FieldTimeSeries` calculated from the salinity and temperature
`FieldTimeSeries` DNS simulation output. The keyword argument `reference_pressure` can be
passed to specify a reference pressure at which to compute the density variable.
"""
function compute_density(S_timeseries::FieldTimeSeries, T_timeseries::FieldTimeSeries;
                         reference_pressure = 0)

    ρ_ts = deepcopy(S_timeseries)
    t = S_timeseries.times
    for i ∈ eachindex(t)
        Sᵢ, Θᵢ = S_timeseries[i], T_timeseries[i]
        ρ_ts[i] .= @at (Center, Center, Center) gsw_rho.(Sᵢ, Θᵢ, reference_pressure)
    end

    return ρ_ts

end

end # module
