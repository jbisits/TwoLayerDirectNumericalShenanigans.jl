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
"""
    function TLDNS.animate_tracers(tracers::AbstractString)
Animate the salinity and temperature `tracers` from saved `.nc` output.
"""
function TLDNS.animate_tracers(tracers::AbstractString; xslice = 52, yslice = 52)

    NCDataset(tracers) do ds

        x = ds["xC"][:]
        z = ds["zC"][:]
        t = ds["time"][:]

        pred_Θₗ, pred_Sₗ = ds.attrib["Predicted equilibrium Tₗ"],
                            ds.attrib["Predicted equilibrium Sₗ"]

        n = Observable(1)
        S = @lift ds["S"][:, yslice, :, $n]
        S_profile = @lift ds["S"][xslice, yslice, :, $n]
        Θ = @lift ds["T"][:, yslice, :, $n]
        Θ_profile = @lift ds["T"][xslice, yslice, :, $n]
        time_title = @lift @sprintf("t=%1.2f minutes", t[$n] / 60)

        fig = Figure(size = (1000, 1000))
        ax = [Axis(fig[j, i], title = (i == 1 && j == 1) ? time_title : "") for i ∈ 1:2, j ∈ 1:2]

        # Salinity
        lines!(ax[1], S_profile, z)
        ax[1].xlabel = "S gkg⁻¹"
        ax[1].ylabel = "z (m)"
        ax[1].xaxisposition = :top
        vlines!(ax[1], pred_Sₗ, linestyle = :dash, color = :red,
                label = "Predicted Sₗ")
        axislegend(ax[1], position = :lb)

        Scmap = cgrad(:haline)[2:end-1]
        Srange = extrema(ds[:S][:, :, :, 1])
        Slow = cgrad(:haline)[1]
        Shigh = cgrad(:haline)[end]
        hm = heatmap!(ax[2], x, z, S, colorrange = Srange, colormap = Scmap,
                        lowclip = Slow, highclip = Shigh)

        ax[2].xlabel = "x (m)"
        ax[2].ylabel = "z (m)"
        Colorbar(fig[1, 3], hm, label = "S gkg⁻¹")

        linkyaxes!(ax[1], ax[2])
        hideydecorations!(ax[2], ticks = false)

        # Temperature
        lines!(ax[3], Θ_profile, z)
        ax[3].xlabel = "Θ°C"
        ax[3].ylabel = "z (m)"
        ax[3].xaxisposition = :top
        vlines!(ax[3], pred_Θₗ, linestyle = :dash, color = :red,
                label = "Predicted Θₗ")
        axislegend(ax[3], position = :lb)

        Θcmap = cgrad(:thermal)[2:end-1]
        Θrange = extrema(ds[:T][:, :, :, 1])
        Θlow = cgrad(:thermal)[1]
        Θhigh = cgrad(:thermal)[end]
        hm = heatmap!(ax[4], x, z, Θ, colorrange = Θrange, colormap = Θcmap,
                        lowclip = Θlow, highclip = Θhigh)

        ax[4].xlabel = "x (m)"
        ax[4].ylabel = "z (m)"
        Colorbar(fig[2, 3], hm, label = "Θ°C")

        linkyaxes!(ax[3], ax[4])
        hideydecorations!(ax[4], ticks = false)

        frames = eachindex(t)
        record(fig, joinpath(pwd(), "tracers.mp4"),
            frames, framerate=8) do i
            msg = string("Plotting frame ", i, " of ", frames[end])
            print(msg * " \r")
            n[] = i
        end

    end

    return nothing
end
"""
    function animate_density(computed_output::AbstractString, variable::AbstractString;
                                     xslice = 52, yslice = 52)
Animate the density `variable` in `computed_output`.
"""
function TLDNS.animate_density(computed_output::AbstractString, variable::AbstractString;
                               xslice = 52, yslice = 52)

    NCDataset(computed_output) do ds

        x = ds["xC"][:]
        z = ds["zC"][:]
        t = ds["time"][:]

        pred_max_density = ds.attrib["Predicted maximum density"]

        n = Observable(1)
        σ = @lift ds[variable][:, yslice, :, $n]
        σ_profile = @lift ds[variable][xslice, yslice, :, $n]
        time_title = @lift @sprintf("t=%1.2f minutes", t[$n] / 60)

        fig = Figure(size = (1000, 500))
        ax = [Axis(fig[1, i], title = i == 1 ? time_title : "") for i ∈ 1:2]

        lines!(ax[1], σ_profile, z)
        ax[1].xlabel = "S gkg⁻¹"
        ax[1].ylabel = "z"
        ax[1].xaxisposition = :top
        vlines!(ax[1], pred_max_density, linestyle = :dash, color = :red,
                label = "Predicted Sₗ")
        axislegend(ax[1], position = :lb)

        colormap = cgrad(:dense)[2:end-1]
        colorrange = (minimum(ds[variable][:, :, :, 1]), pred_max_density)
        lowclip = cgrad(:dense)[1]
        highclip = cgrad(:dense)[end]
        hm = heatmap!(ax[2], x, z, σ; colorrange, colormap, lowclip, highclip)

        ax[2].xlabel = "x (m)"
        ax[2].ylabel = "z (m)"
        Colorbar(fig[1, 3], hm, label = "σ₀ kgm⁻³")

        linkyaxes!(ax[1], ax[2])
        hideydecorations!(ax[2], ticks = false)

        frames = eachindex(t)
        record(fig, joinpath(pwd(), "density.mp4"),
            frames, framerate=8) do i
            msg = string("Plotting frame ", i, " of ", frames[end])
            print(msg * " \r")
            n[] = i
        end

    end

    return nothing
end
"""
    function plot_scalar_diagnostics(computed_output::AbstractString)
Animate the density and diagnostics (at this stage ∫ϵ and ∫κᵥ) from `computed_output`. **Note:**
this function assumes that `computed_output` is a `.nc` file.
"""
function TLDNS.plot_scalar_diagnostics(computed_output::AbstractString)

    NCDataset(computed_output) do ds
        t = ds["time"][:] / (60)
        ∫ϵ = ds["∫ϵ"][:]
        ∫κᵥ = ds["∫κᵥ"][:]

        fig = Figure(size = (600, 1000))
        ax = [Axis(fig[i, 1]) for i ∈ 1:2]
        lines!(ax[1], t, ∫ϵ)
        ax[1].title = "Integraged TKE"
        ax[1].ylabel = "∫ϵ"

        lines!(ax[2], t, ∫κᵥ)
        ax[2].title = "Integrated inferred vertical diffusivity"
        ax[2].xlabel = "t (minutes)"
        ax[2].ylabel = "∫κᵥ"

        linkxaxes!(ax[1], ax[2])

        plotsave = "TKE_and_inferredK.png"
        save(plotsave, fig)
        @info "Save plot to $(plotsave)"
    end

    return nothing

end
"""
    function hovmoller(computed_output::AbstractString, variable::AbstractString)
Produce a Hovmoller plot of a `variable` in `computed_output`.
"""
function TLDNS.hovmoller(computed_output::AbstractString, variable::AbstractString;
                         colormap = :viridis, unit = nothing)

    NCDataset(computed_output) do ds

        t = ds["time"][:] ./ 60
        zC = ds["zC"][:]
        plot_var = ds[variable][:, :]
        fig = Figure(size = (600, 1000))
        ax = Axis(fig[1, 1],
                  title = "Hovmoller plot of $(variable)",
                  xlabel = "t (minutes)",
                  ylabel = "z (m)")
        colormap = cgrad(colormap)[2:end-1]
        lowclip = cgrad(colormap)[1]
        highclip = cgrad(colormap)[end]
        hm = heatmap!(ax, t, zC, plot_var';
                      colormap, lowclip, highclip, colorrange = (-1e-4, 1e-4))
        cbar_label = isnothing(unit) ? variable : variable * unit
        Colorbar(fig[1, 2], hm, label = cbar_label)

        plotsave = "hovmoller_"*variable*".png"
        save(plotsave, fig)
        @info "Save plot to $(plotsave)"

    end

    return nothing

end
"""
    function animate_tracer_distributions(tracers::AbstractString;
                                          S_binwidth = 0.001, Θ_binwidth = 0.01)
Animate volume distribution of the salinity and temperature in `tracers`.
**Note:** this assumens `tracers` is a `.nc` file.
"""
function TLDNS.animate_tracer_distributions(tracers::AbstractString;
                                            S_binwidth = 0.001, Θ_binwidth = 0.01)

    NCDataset(tracers) do ds

        t = ds["time"][:]

        pred_Θₗ, pred_Sₗ = ds.attrib["Predicted equilibrium Tₗ"],
                            ds.attrib["Predicted equilibrium Sₗ"]

        n = Observable(1)
        S_extrema = extrema(ds["S"][:, :, :, 1])
        S_edges = S_extrema[1]-S_binwidth:S_binwidth:S_extrema[2]+S_binwidth
        S_hist = @lift fit(Histogram, reshape(ds["S"][:, :, :, $n], :), S_edges)
        Θ_extrema = extrema(ds["T"][:, :, :, 1])
        Θ_edges = Θ_extrema[1]-Θ_binwidth:Θ_binwidth:Θ_extrema[2]+Θ_binwidth
        Θ_hist = @lift fit(Histogram, reshape(ds["T"][:, :, :, $n], :), Θ_edges)
        time_title = @lift @sprintf("t=%1.2f minutes", t[$n] / 60)

        fig = Figure(size = (1000, 500))
        ax = [Axis(fig[1, i], title = i == 1 ? time_title : "") for i ∈ 1:2]

        hist!(ax[1], S_hist)
        ax[1].xlabel = "S (gkg⁻¹)"
        ax[1].ylabel = "Frequency"
        vlines!(ax[1], pred_Sₗ)
        hist!(ax[2], Θ_hist)
        ax[2].xlabel = "Θ (°C)"
        ax[2].ylabel = "Frequency"
        vlines!(ax[2], pred_Θₗ)

        hideydecorations!(ax[2], ticks = false)
        linkyaxes!(ax[1], ax[2])

        frames = eachindex(t)
        record(fig, joinpath(pwd(), "S_and_T_distributions.mp4"),
            frames, framerate=8) do i
            msg = string("Plotting frame ", i, " of ", frames[end])
            print(msg * " \r")
            n[] = i
        end

    end

    return nothing

end
"""
    function animate_joint_tracer_distribution(tracers::AbstractString;
                                               S_binwidth = 0.001, Θ_binwidth = 0.01)
Animate the joint distribution of the salinity and temperature in `tracers`.
**Note:** this assumens `tracers` is a `.nc` file.
"""
function TLDNS.animate_joint_tracer_distribution(tracers::AbstractString;
                                                 S_binwidth = 0.001, Θ_binwidth = 0.01)

    NCDataset(tracers) do ds

        t = ds["time"][:]

        pred_Θₗ, pred_Sₗ = ds.attrib["Predicted equilibrium Tₗ"],
                            ds.attrib["Predicted equilibrium Sₗ"]

        n = Observable(1)
        S_extrema = extrema(ds["S"][:, :, :, 1])
        S_edges = S_extrema[1]-S_binwidth:S_binwidth:S_extrema[2]+S_binwidth
        Θ_extrema = extrema(ds["T"][:, :, :, 1])
        Θ_edges = Θ_extrema[1]-Θ_binwidth:Θ_binwidth:Θ_extrema[2]+Θ_binwidth
        joint_hist = @lift fit(Histogram, (reshape(ds["S"][:, :, :, $n], :),
                                           reshape(ds["T"][:, :, :, $n], :)),
                                          (S_edges, Θ_edges))
        time_title = @lift @sprintf("t=%1.2f minutes", t[$n] / 60)

        fig = Figure(size = (500, 500))
        ax = Axis(fig[1, 1], title =  time_title)

        hm = heatmap!(ax, S_edges, Θ_edges, replace(joint_hist.weights, 0 => NaN))
        scatter!(ax, [pred_Sₗ], [pred_Θₗ], color = :red, label = "Predcited (Sₗ, Θₗ)")
        ax.xlabel = "S (gkg⁻¹)"
        ax.ylabel = "Θ (°C)"
        Colorbar(fig[1, 2], hm)
        axislegend(ax, position = :lt)

        frames = eachindex(t)
        record(fig, joinpath(pwd(), "S_T_joint_distribution.mp4"),
        frames, framerate=8) do i
            msg = string("Plotting frame ", i, " of ", frames[end])
            print(msg * " \r")
            n[] = i
        end

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
