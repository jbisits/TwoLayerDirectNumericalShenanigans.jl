"""
    function animate_2D_field(field_timeseries::FieldTimeSeries, field_name::AbstractString,
                              field_dimensions::NamedTuple{Symbol, Symbol}; colormap = :thermal)
Animate a time series that is saved in memory.

Function arguments:

- `field_timeseries` to be animated;
- `field_name` the name of the variable in the `field_timeseries`;
- the `xslice` and `yslice` for the `x-z` heatmap and vertical profile.

Keyword arguments:

- `colormap` for the animated `field_timeseries`;
- `aspect_ration` for the animation;
- `z_extrema` where to look for the extreme z values while avoiding any random noise, pass
two ranges or two values for where to look for the `extrema`.
"""
function TLDNS.animate_2D_field(field_timeseries::FieldTimeSeries, field_name::AbstractString,
                               xslice::Int64, yslice::Int64; colormap = :thermal,
                               colorrange = nothing, highclip = nothing, lowclip = nothing,
                               aspect_ratio = 1)

    x, y, z = nodes(field_timeseries[1])

    t = field_timeseries.times
    n = Observable(1)
    field_tₙ = @lift interior(field_timeseries[$n], :, yslice, :)
    profile_tₙ = @lift interior(field_timeseries[$n], xslice, yslice, :)
    colorrange = isnothing(colorrange) ? extrema(interior(field_timeseries, :, :, :, 1)) :
                                         colorrange
    time_title = @lift @sprintf("t=%1.2f minutes", t[$n] / 60)

    fig = Figure(size = (1000, 600))
    ax = [Axis(fig[1, i], title = i == 1 ? time_title : "") for i ∈ 1:2]

    hm = heatmap!(ax[1], x, z, field_tₙ; colorrange, colormap, lowclip, highclip)

    ax[1].xlabel = "x (m)"
    ax[1].ylabel = "z (m)"
    ax[1].aspect = aspect_ratio
    ax[1].xticklabelrotation = π / 4
    ax[1].xlabel = "x"
    ax[1].ylabel = "y"
    Colorbar(fig[2, 1], hm, label = field_name, vertical = false, flipaxis = false)

    lines!(ax[2], profile_tₙ, z)
    ax[2].xlabel = field_name
    ax[2].ylabel = "z"
    ax[2].aspect = aspect_ratio
    ax[2].xaxisposition = :top

    linkyaxes!(ax[1], ax[2])

    frames = eachindex(t)
    record(fig, joinpath(pwd(), field_name * ".mp4"),
          frames, framerate=8) do i
        msg = string("Plotting frame ", i, " of ", frames[end])
        print(msg * " \r")
        n[] = i
    end

    return nothing

end
function TLDNS.animate_2D_field(rs::Raster, xslice::Int64, yslice::Int64; colormap = :thermal,
                                colorrange = nothing, highclip = nothing, lowclip = nothing,
                                aspect_ratio = 1)

    x, z, t = lookup(rs, :xC), lookup(rs, :zC), lookup(rs, :Ti)
    field_name = string(rs.name)

    n = Observable(1)
    field_tₙ = @lift rs.data[:, yslice, :, $n]
    profile_tₙ = @lift rs.data[xslice, yslice, :, $n]
    colorrange = isnothing(colorrange) ? extrema(rs) : colorrange
    time_title = @lift @sprintf("t=%1.2f minutes", t[$n] / 60)

    fig = Figure(size = (1000, 600))
    ax = [Axis(fig[1, i], title = i == 1 ? time_title : "") for i ∈ 1:2]

    hm = heatmap!(ax[1], x, z, field_tₙ; colorrange, colormap, lowclip, highclip)

    ax[1].xlabel = "x (m)"
    ax[1].ylabel = "z (m)"
    ax[1].aspect = aspect_ratio
    ax[1].xticklabelrotation = π / 4
    ax[1].xlabel = "x"
    ax[1].ylabel = "y"
    Colorbar(fig[2, 1], hm, label = field_name, vertical = false, flipaxis = false)

    lines!(ax[2], profile_tₙ, z)
    ax[2].xlabel = field_name
    ax[2].ylabel = "z"
    ax[2].aspect = aspect_ratio
    ax[2].xaxisposition = :top

    linkyaxes!(ax[1], ax[2])

    frames = eachindex(t)
    record(fig, joinpath(pwd(), field_name * ".mp4"),
        frames, framerate=8) do i
        msg = string("Plotting frame ", i, " of ", frames[end])
        print(msg * " \r")
        n[] = i
    end

    return nothing

end
"""
    function visualise_initial_stepchange(dns::TwoLayerDNS, interface_location::Number)
Plot an initial step change of the `tracers` in a `model`. This function assumes there are two
tracers (salinity and temperature) and plots the x-z, y-z and field-z initial fields.
"""
function TLDNS.visualise_initial_stepchange(dns::TwoLayerDNS, interface_location::Number)

    model, initial_conditions = dns.model, dns.initial_conditions
    x, y, z = nodes(model.grid, (Center(), Center(), Center()))
    S₀ˡ, ΔS₀ = initial_conditions.S₀ˡ, initial_conditions.ΔS₀
    T₀ˡ, ΔT₀ = initial_conditions.T₀ˡ, initial_conditions.ΔT₀
    S_stepchange = initial_tracer_heaviside.(z, S₀ˡ, ΔS₀, interface_location)
    S_sc_array = repeat(S_stepchange', length(x), 1)
    T_stepchange = initial_tracer_heaviside.(z, T₀ˡ, ΔT₀, interface_location)
    T_sc_array = repeat(T_stepchange', length(x), 1)
    fig = Figure(size = (600, 1500))
    ax = [Axis(fig[j, i]) for i ∈ 1:2, j ∈ 1:3]

    hm = heatmap!(ax[1], x, z, S_sc_array; colormap = :haline)
    ax[1].title = "Initial salinity (x-z)"
    ax[1].xlabel = "x (m)"
    ax[1].ylabel = "z (m)"
    heatmap!(ax[2], x, z, S_sc_array; colormap = :haline)
    ax[2].title = "Initial salinity (y-z)"
    ax[2].xlabel = "y (m)"
    ax[2].ylabel = "z (m)"
    Colorbar(fig[1, 3], hm, label = "S (gkg⁻¹)")
    hm = heatmap!(ax[3], y, z, T_sc_array; colormap = :thermal)
    ax[3].title = "Initial temperature (x-z)"
    ax[3].xlabel = "x (m)"
    ax[3].ylabel = "z (m)"
    heatmap!(ax[4], y, z, T_sc_array; colormap = :thermal)
    ax[4].title = "Initial temperature (y-z)"
    ax[4].xlabel = "y (m)"
    ax[4].ylabel = "z (m)"
    Colorbar(fig[2, 3], hm, label = "Θ (°C)")
    lines!(ax[5], S_stepchange, z)
    ax[5].title = "Initial salinity profile"
    ax[5].xlabel = "S (gkg⁻¹)"
    ax[5].ylabel = "z (m)"
    lines!(ax[6], T_stepchange, z)
    ax[6].title = "Initial temperature profile"
    ax[6].xlabel =  "Θ (°C)"
    ax[6].ylabel = "z (m)"

    return fig

end
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
TLDNS.initial_tracer_heaviside(z, C::Number, ΔC::Number, interface_location) =
                                                    z - interface_location < 0 ? C : C + ΔC
"""
    function visualise_initial_conditions(dns::TwoLayerDNS, xslice::Integer, yslice::Integer)
Plot the initial state of the `tracers` in a `model`. This function assumes there are two
tracers (salinity and temperature) and plots the x-z, y-z and field-z initial fields at
`xslice` and `yslice`.
"""
function TLDNS.visualise_initial_conditions(dns::TwoLayerDNS, xslice::Integer, yslice::Integer)

    model = dns.model

    x, y, z = nodes(model.grid, (Center(), Center(), Center()))
    S = model.tracers.S
    T = model.tracers.T
    fig = Figure(size = (600, 1500))
    ax = [Axis(fig[j, i]) for i ∈ 1:2, j ∈ 1:3]

    hm = heatmap!(ax[1], x, z, interior(S, :, yslice, :, 1); colormap = :haline)
    ax[1].title = "Initial salinity (x-z)"
    ax[1].xlabel = "x (m)"
    ax[1].ylabel = "z (m)"
    heatmap!(ax[2], x, z, interior(S, xslice, :, :, 1); colormap = :haline)
    ax[2].title = "Initial salinity (y-z)"
    ax[2].xlabel = "y (m)"
    ax[2].ylabel = "z (m)"
    Colorbar(fig[1, 3], hm, label = "S (gkg⁻¹)")
    hm = heatmap!(ax[3], x, z, interior(T, :, yslice, :, 1); colormap = :thermal)
    ax[3].title = "Initial temperature (x-z)"
    ax[3].xlabel = "x (m)"
    ax[3].ylabel = "z (m)"
    heatmap!(ax[4], y, z, interior(T, xslice, :, :, 1); colormap = :thermal)
    ax[4].title = "Initial temperature (y-z)"
    ax[4].xlabel = "y (m)"
    ax[4].ylabel = "z (m)"
    Colorbar(fig[2, 3], hm, label = "Θ (°C)")
    lines!(ax[5], interior(S, xslice, yslice, :, 1), z)
    ax[5].title = "Initial salinity profile"
    ax[5].xlabel = "S (gkg⁻¹)"
    ax[5].ylabel = "z (m)"
    lines!(ax[6], interior(T, xslice, yslice, :, 1), z)
    ax[6].title = "Initial temperature profile"
    ax[6].xlabel =  "Θ (°C)"
    ax[6].ylabel = "z (m)"

    return fig

end
"""
    function visualise_initial_density(dns::TwoLayerDNS, xslice::Integer,  yslice::Integer,
                                       pressure::Union{Number, Vector{Number}})
Compute and plot the initial density at `pressure` (either reference pressure or in-situ
pressure). The arguments `xslice` and `yslice` are used to choose where in the domain the
figures are from.
"""
function TLDNS.visualise_initial_density(dns::TwoLayerDNS, xslice::Integer,  yslice::Integer,
                                        pressure::Union{Number, Vector{Number}})

    model = dns.model

    x = xnodes(model.grid, Center(), Center(), Center())
    z = znodes(model.grid, Center(), Center(), Center())
    S = interior(model.tracers.S, :, yslice, :, 1)
    T = interior(model.tracers.T, :, yslice, :, 1)
    ρ = gsw_rho.(S, T, pressure)

    fig = Figure(size = (1000, 600))
    ax = [Axis(fig[1, i]) for i ∈ 1:2]

    hm = heatmap!(ax[1], x, z, ρ; colormap = :dense)
    ax[1].title = "Initial density (x-z)"
    ax[1].xlabel = "x (m)"
    ax[1].ylabel = "z (m)"
    Colorbar(fig[2, 1], hm, label = "ρ (kgm⁻³)", vertical = false, flipaxis = false)
    lines!(ax[2], ρ[xslice, :], z)
    ax[2].title = "Initial density profile"
    ax[2].xlabel = "ρ (kgm⁻³)"
    ax[2].ylabel = "z (m)"

    linkyaxes!(ax[1], ax[2])

    return fig

end

"""
function visualise_snapshot(field_timeseries::FieldTimeSeries, field_name::AbstractString,
                            snapshot::Int64)
Plot a `snapshot` of the `field_timeseries`  with `field_name` at `xslice`, `yslice`.
"""
function TLDNS.visualise_snapshot(field_timeseries::FieldTimeSeries, field_name::AbstractString,
                                 xslice::Int64, yslice::Int64, snapshot::Int64;
                                 colormap = :thermal)

    x, y, z = nodes(field_timeseries[1])
    t = round(field_timeseries.times[snapshot]; digits = 3)
    fig = Figure(size = (1000, 500))
    ax = [Axis(fig[1, i]) for i ∈ 1:2]

    hm = heatmap!(ax[1], x, z, interior(field_timeseries, :, yslice, :, snapshot); colormap)
    ax[1].title = field_name * " at time t = $(t) (x-z)"
    ax[1].xlabel = "x (m)"
    ax[1].ylabel = "z (m)"
    Colorbar(fig[2, 1], hm, vertical = false, label = field_name, flipaxis = false)
    lines!(ax[2], interior(field_timeseries, xslice, yslice, :, snapshot), z)
    ax[2].title = field_name * " profile at time t = $(t)"
    ax[2].xlabel = field_name
    ax[2].ylabel = "z (m)"

    linkyaxes!(ax[1], ax[2])

    return fig

end

function TLDNS.visualise_snapshot(rs::Raster, yslice::Int64, snapshot::Int64;
                                  colormap = :thermal, unit = nothing, aspect_ratio = 1)

        x, z, t = lookup(rs, :xC), lookup(rs, :zC), lookup(rs, :Ti)
        field_name = string(rs.name)
        fig = Figure(size = (500, 500))
        ax = Axis(fig[1, 1])

        hm = heatmap!(ax, x, z, rs.data[:, yslice, :, snapshot]; colormap)
        ax.title = field_name * " at time t = $(t) (x-z)"
        ax.xlabel = "x (m)"
        ax.ylabel = "z (m)"
        ax.aspect = aspect_ratio
        cbar_label = isnothing(unit) ? field_name : field_name * unit
        Colorbar(fig[1, 2], hm, label = cbar_label)

    return fig

end
function TLDNS.visualise_snapshot(rs::Raster, xslice::Int64, yslice::Int64, snapshot::Int64;
                                  colormap = :thermal, aspect_ratio = 1)

    x, z, t = lookup(rs, :xC), lookup(rs, :zC), lookup(rs, :Ti)
    field_name = string(rs.name)
    fig = Figure(size = (1000, 500))
    ax = [Axis(fig[1, i]) for i ∈ 1:2]

    hm = heatmap!(ax[1], x, z, rs.data[:, yslice, :, snapshot]; colormap)
    ax[1].title = field_name * " at time t = $(t) (x-z)"
    ax[1].xlabel = "x (m)"
    ax[1].ylabel = "z (m)"
    ax[1].aspect = aspect_ratio
    Colorbar(fig[2, 1], hm, vertical = false, label = field_name, flipaxis = false)
    lines!(ax[2], rs.data[xslice, yslice, :, snapshot], z)
    ax[2].title = field_name * " profile at time t = $(t)"
    ax[2].xlabel = field_name
    ax[2].ylabel = "z (m)"

    linkyaxes!(ax[1], ax[2])

    return fig

end
