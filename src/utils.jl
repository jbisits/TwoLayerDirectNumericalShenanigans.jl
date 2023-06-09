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
