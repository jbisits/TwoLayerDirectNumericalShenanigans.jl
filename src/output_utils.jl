"""
    function compute_density(S_timeseries::FieldTimeSeries, T_timeseries::FieldTimeSeries;
                             reference_pressure = 0)
Return a density `FieldTimeSeries` calculated from the salinity and temperature
`FieldTimeSeries` DNS simulation output. The keyword argument `reference_pressure` can be
passed to specify a reference pressure at which to compute the density variable.
"""
function compute_density(S_timeseries::FieldTimeSeries, T_timeseries::FieldTimeSeries;
                         reference_pressure = 0)

    ρ_ts = FieldTimeSeries{Center, Center, Center}(S_timeseries.grid, S_timeseries.times,
                                                  indices = S_timeseries.indices,
                                                  boundary_conditions =
                                                    S_timeseries.boundary_conditions)

    t = S_timeseries.times
    for i ∈ eachindex(t)
        Sᵢ, Θᵢ = S_timeseries[i], T_timeseries[i]
        ρ_ts[i] .= @at (Center, Center, Center) gsw_rho.(Sᵢ, Θᵢ, reference_pressure)
    end

    return ρ_ts

end
"""
    function compute_density!(filepath::String; reference_pressure = 0)
Compute a density variable at `reference_pressure` and append it to saved output. **Note**
this only needs to be done once! After that the group will exist and will not be able to be
overwritten but a different group can be made by passing a different `density_string`.
This allows calling `FieldTimeSeries` on the `density_variable`
"""
function compute_density!(filepath::AbstractString; density_string = "σ", reference_pressure = 0)

    file_type = find_file_type(filepath)
    if isequal(file_type, ".nc")

        @info "Lazily loading S and T into separate Rasters"
        S_rs = Raster(filepath, lazy = true, name = :S)
        T_rs = Raster(filepath, lazy = true, name = :T)
        time = lookup(S_rs, :Ti)

        NCDataset(filepath, "a") do ds
            @info "Appending density variable to saved .nc file"
            defVar(ds, density_string, zeros(size(S_rs)), ("xC", "yC", "zC", "time"),
                    attrib = Dict("units" => "kgm⁻³",
                                  "longname" => "Potential density",
                                  "comments" => "computed at reference pressues p = $reference_pressure"))
            @info "Computing density and saving to .nc file"
            for t ∈ eachindex(time)
                ds[density_string][:, :, :, t] = get_σₚ(S_rs[:, :, :, t], T_rs[:, :, :, t],
                                                        reference_pressure)
            end
        end

    elseif isequal(file_type, ".jld2")

        file = jldopen(filepath, "a+")
        file_keys = keys(file["timeseries"]["S"])
        JLD2.Group(file, "timeseries/"*density_string)
        for (i, key) ∈ enumerate(file_keys)
            if i == 1
                for k ∈ keys(file["timeseries/S/"*key])
                    file["timeseries/"*density_string*"/"*key*"/"*k] =
                        file["timeseries/S/"*key*"/"*k]
                end
                file["timeseries/"*density_string*"/"*key*"/reference_pressure"] =
                    reference_pressure
            else
                Sᵢ, Θᵢ = file["timeseries/S/"*key], file["timeseries/T/"*key]
                file["timeseries/"*density_string*"/"*key]= gsw_rho.(Sᵢ, Θᵢ, reference_pressure)
            end
        end

        close(file)

    end

    return nothing

end
"""
    function append_density!(; saved_simulations = readdir(SIMULATION_PATH, join = true))
Append a density timeseries to the saved output from a `TwoLayerDNS` simulation. By default
the function looks for saved data at the default filepath for saving `SIMULATION_PATH`.
Pass another path as the keyword argument `saved_simulations` to look somewhere else.
"""
function append_density!(; saved_simulations = readdir(SIMULATION_PATH, join = true))

    for simulation ∈ saved_simulations

        file_type = find_file_type(simulation)

        if isequal(file_type, ".jld2")
            open_sim = jldopen(simulation)
            if "σ" ∉ keys(open_sim["timeseries"])
                close(open_sim)
                compute_density!(simulation, density_string = "σ₀")
            else
                @info "A density timeseries already exists in $simulation."
                close(open_sim)
            end
        elseif isequal(file_type, ".nc")
            NCDataset(simulation, "a") do ds
                if "σ" ∉ keys(ds)
                    compute_density!(simulation)
                else
                    @info "A density timeseries already exists in $simulation."
                end
            end
        end
    end

    return nothing

end

"Return the mean from a `FieldTimeSeries` that is `OnDisk()`."
function field_ts_timemean(field_ts::FieldTimeSeries)

    t = field_ts.times
    field_data = field_ts[1].data
    for i ∈ 2:length(t)
        field_data .+= field_ts[i].data
    end

    return field_data ./ length(t)

end
"""
    function predicted_maximum_density()
Compute the predicted maximum density of water that can form along the mixing line.
"""
function predicted_maximum_density!(simulation::Simulation, dns::TwoLayerDNS; reference_pressure = 0)

    Sᵘ, Sˡ, ΔS, Tᵘ, Tˡ, ΔT = dns.initial_conditions
    slope = - ΔT / ΔS
    S_mix = range(Sᵘ, Sˡ, step = 0.000001)
    T_mix = @. Tᵘ + slope * (Sᵘ - S_mix)
    ρ_mix_max = maximum(gsw_rho.(S_mix, T_mix, reference_pressure))
    NCDataset(simulation.output_writers[:outputs].filepath, "a") do ds
        ds.attrib["Predicted maximum density"] = ρ_mix_max
    end

    return nothing

end
function predicted_maximum_density!(file::AbstractString; reference_pressure = 0)

    NCDataset(file, "a") do ds
        Sᵘ, Sˡ = ds["S"][end, end, end, 1], ds["S"][1, 1, 1, 1]
        ΔS = Sᵘ - Sˡ
        Tᵘ, Tˡ = ds["T"][end, end, end, 1], ds["T"][1, 1, 1, 1]
        ΔT = Tᵘ - Tˡ
        slope = - ΔT / ΔS
        S_mix = range(Sᵘ, Sˡ, step = 0.000001)
        T_mix = @. Tᵘ + slope * (Sᵘ - S_mix)
        ρ_mix_max = maximum(gsw_rho.(S_mix, T_mix, reference_pressure))
        ds.attrib["Predicted maximum density"] = ρ_mix_max
    end

    return nothing

end
"Calculate the Kolmogorov length scale `η` from viscousity and average TKE dissapation."
η(ν, ϵ) = (ν^3 / ϵ)^(1/4)
"""
    function minimum_η(ϵ::FieldTimeSeries; ν = 1e-6)
Find the minimum `η`, i.e. the Kolmogorov length scale, from the `KineticEnergyDissaption`,
 `ϵ`, time series.
"""
function minimum_η(ϵ::FieldTimeSeries; ν = 1e-6)

    t = ϵ.times
    minimum_η_t = similar(t)
    for i ∈ eachindex(t)
        minimum_η_t[i] = minimum(η.(ν, ϵ[i].data))
    end

    return minimum(minimum_η_t)

end
function minimum_η(ϵ::Raster; ν = 1e-6)

    t = lookup(ϵ, :Ti)
    minimum_η_t = similar(t)
    for i ∈ eachindex(t)
        minimum_η_t[i] = minimum(η.(ν, ϵ.data[:, :, :, i]))
    end

    return minimum(minimum_η_t)

end
"""
    function kolmogorov_and_batchelor_scale!(file::AbstractString)
Append the minimum Kolmogorov and Batchelor scales (in space and time) from a `TwoLayerDNS`
simulation with output saved on `file`. The Kolmogorov scale is defined by
```math
    η = \\left(\\frac{ν³}{ϵ}\\right)^\\frac{1}{4}
```
and the Batchelor scale is
```math
    λ_{B} = \\frac{η}{√Sc}
```
where ``Sc`` is the Schmidt number.
"""
function kolmogorov_and_batchelor_scale!(file::AbstractString)

    file_type = find_file_type(file)
    if isequal(file_type, ".nc")

        ϵ = Raster(file, lazy = true, name = :ϵ)
        ds = NCDataset(file, "a")
        ν_str = ds.attrib["ν"]
        ν = parse(Float64, ν_str[1:findfirst('m', ν_str)-1])
        Sc = ds.attrib["Sc"]
        η_min = minimum_η(ϵ; ν)
        ds.attrib["η (min)"] = η_min # minimum space and time Kolmogorov scale
        ds.attrib["λ_B"] = η_min / sqrt(Sc) # minimum space and time Batchelor scale
        close(ds)

    elseif isequal(file_type, ".jld2")

        ϵ_ts = FieldTimeSeries(file, "ϵ", backend = OnDisk())
        Sc = load(file, "Non_dimensional_numbers")["Sc"]
        ν = load(file, "closure/ν")
        min_η = minimum_η(ϵ_ts; ν)

        jldopen(file, "a+") do f
            f["minimum_kolmogorov_scale"] = min_η
            f["minimum_batchelor_scale"] = min_η / sqrt(Sc)
        end

    end

    return nothing

end
"""
    function non_dimensional_numbers(simulation::Simulation, dns::TwoLayerDNS)
Compute and append non-dimensional numbers related to the DNS experiments.
The non-dimensional numbers are:

- Prandtl number: ``Pr = ν / κₜ``
- Schmidt number: ``Sc = ν / κₛ``
- Lewis number:   ``Le = κₜ / κₛ``
- Raleigh number (density): ``Ra_{d} = Ra_{t} / Ra_{s} = (αΔT / βΔS) * (1 / Le)``.

These numbers are then saved into the simulation output file.
"""
function non_dimensional_numbers!(simulation::Simulation, dns::TwoLayerDNS)

    model, initial_conditions = dns.model, dns.initial_conditions
    ν = model.closure.ν
    κₛ, κₜ = model.closure.κ
    Pr = ν / κₜ
    Sc = ν / κₛ
    Le = κₜ / κₛ
    α = gsw_alpha(initial_conditions.S₀ˡ, initial_conditions.T₀ˡ, 0)
    β = gsw_beta(initial_conditions.S₀ˡ, initial_conditions.T₀ˡ, 0)
    Ra = ((α * initial_conditions.ΔT₀ )/ (β * initial_conditions.ΔS₀)) * (1 / Le)

    nd_nums = Dict("Pr" => Pr, "Sc" => Sc, "Le" => Le, "Ra_ρ" => Ra)

    if simulation.output_writers[:outputs] isa NetCDFOutputWriter

        ds = NCDataset(simulation.output_writers[:outputs].filepath, "a")
        ds.attrib["EOS"] = summary(model.buoyancy.model.equation_of_state.seawater_polynomial)
        ds.attrib["Reference density"] = "$(model.buoyancy.model.equation_of_state.reference_density)kgm⁻³"
        ds.attrib["ν"]  = "$(model.closure.ν) m²s⁻¹"
        ds.attrib["κₛ"] = "$(model.closure.κ.S) m²s⁻¹"
        ds.attrib["κₜ"] = "$(model.closure.κ.T) m²s⁻¹"
        for key ∈ keys(nd_nums)
            ds.attrib[key] = nd_nums[key]
        end
        close(ds)

    else

        jldopen(simulation.output_writers[:outputs].filepath, "a+") do file
            file["Non_dimensional_numbers"] = nd_nums
        end

    end

    return nothing

end
"""
    funciton find_file_type(file::AbstractString)
Return the file type (either `.nc` or `.jld2`) of a `file`.
"""
find_file_type(file::AbstractString) = file[findlast('.', file):end]
