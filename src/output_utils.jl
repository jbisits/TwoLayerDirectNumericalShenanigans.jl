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
function compute_density!(filepath::String; reference_pressure = 0, density_string = "ρ")

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
        open_sim = jldopen(simulation)
        if "σ₀" ∉ keys(open_sim["timeseries"])
            close(open_sim)
            compute_density!(simulation, density_string = "σ₀")
        else
            @info "A density timeseries already exists in $simulation."
            close(open_sim)
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
