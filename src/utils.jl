"""
    function animate_2D_field(field_timeseries::FieldTimeSeries, field_name::AbstractString,
                              field_dimensions::NamedTuple, file_name::AbstractString;
                              colormap = :thermal)
Animate a time series that is saved in memory. Pass in the timeseries, the name of the field,
a `NamedTuple` of the dimensions along and optional colourmap. Filename is inferred from the
`Field` and the dimensions.
"""
function animate_2D_field(field_timeseries::FieldTimeSeries, field_name::AbstractString,
                          field_dimensions::NamedTuple; colormap = :thermal)

    t = field_timeseries.times
    n = Observable(1)
    field_tₙ = @lift interior(field_timeseries[$n], :, 1, :)
    c_limits = extrema(interior(field_timeseries, :, :, :, 1))
    title = @lift @sprintf("t=%1.2f", t[$n])
    fig, ax, hm = heatmap(field_dimensions[1], field_dimensions[2], field_tₙ;
                          colormap, colorrange = c_limits)
    ax.xlabel = string(keys(field_dimensions)[1])
    ax.ylabel = string(keys(field_dimensions)[2])
    Colorbar(fig[1, 2], hm, label = field_name)

    frames = eachindex(t)
    filename = string(keys(field_dimensions)[1]) * string(keys(field_dimensions)[2]) * "_" *
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
    function compute_density(S_timeseries::FieldTimeSeries, T_timeseries::FieldTimeSeries;
                             reference_pressure = 0)
Return a density `FieldTimeSeries` calculated from the salinity and temperature
`FieldTimeSeries` DNS simulation output. The keyword argument `reference_pressure` can be
passed to specify a reference pressure at which to compute the density variable.
"""
function compute_density(S_timeseries::FieldTimeSeries, T_timeseries::FieldTimeSeries;
                         reference_pressure = 0)

    ρ_ts = deepcopy(S_timeseries)
    for i ∈ eachindex(t)
        Sᵢ, Θᵢ = S_timeseries[i], T_timeseries[i]
        ρ_ts[i] .= @at (Center, Center, Center) gsw_rho.(Sᵢ, Θᵢ, reference_pressure)
    end

    return ρ_ts

end
