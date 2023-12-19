"""
    function predicted_maximum_density
Compute the predicted maximum density of water that can form along the mixing line between
the salinity and temperature of the upper and lower layers.
"""
function predicted_maximum_density!(simulation::Simulation, dns::TwoLayerDNS; reference_gp_height = 0)

    Sᵘ, Sˡ, ΔS, Tᵘ, Tˡ, ΔT = dns.initial_conditions
    slope = ΔT / ΔS
    S_mix = range(Sᵘ, Sˡ, step = 0.000001)
    T_mix = @. Tᵘ + slope * (S_mix - Sᵘ)
    Zᵣ = similar(S_mix)
    fill!(Zᵣ, reference_gp_height)
    eos = dns.model.buoyancy.model.equation_of_state
    ρ_mix_max, max_idx = findmax(ρ.(T_mix, S_mix, Zᵣ, fill(eos, length(S_mix))))
    file_type = find_file_type(simulation.output_writers[:tracers].filepath)
    if isequal(file_type, ".nc")
        NCDataset(simulation.output_writers[:computed_output].filepath, "a") do ds
            ds.attrib["Predicted maximum density"] = ρ_mix_max
        end
        NCDataset(simulation.output_writers[:tracers].filepath, "a") do ds
            ds.attrib["Predicted equilibrium Tₗ"] = T_mix[max_idx]
            ds.attrib["Predicted equilibrium Sₗ"] = S_mix[max_idx]
        end
    elseif isequal(file_type, ".jld2")
        jldopen(simulation.output_writers[:computed_output].filepath, "a+") do f
            f["Predicted maximum density"] = ρ_mix_max
        end
        jldopen(simulation.output_writers[:tracers].filepath, "a+") do f
            f["Predicted equilibrium Tₗ"] = T_mix[max_idx]
            f["Predicted equilibrium Sₗ"] = S_mix[max_idx]
        end
    end

    return nothing

end
function predicted_maximum_density!(file::AbstractString; reference_pressure = 0)

    NCDataset(file, "a") do ds
        Sᵘ, Sˡ = ds["S"][end, end, end, 1], ds["S"][1, 1, 1, 1]
        ΔS = Sᵘ - Sˡ
        Tᵘ, Tˡ = ds["T"][end, end, end, 1], ds["T"][1, 1, 1, 1]
        ΔT = Tᵘ - Tˡ
        slope = ΔT / ΔS
        S_mix = range(Sᵘ, Sˡ, step = 0.000001)
        T_mix = @. Tᵘ + slope * (S_mix - Sᵘ)
        # Do not have eos so cannot use `SeawaterPolynomials.ρ`
        ρ_mix_max, max_idx = findmax(gsw_rho.(S_mix, T_mix, reference_pressure))
        ds.attrib["Predicted maximum density"] = ρ_mix_max
        ds.attrib["Predicted equilibrium Tₗ"] = T_mix[max_idx]
        ds.attrib["Predicted equilibrium Sₗ"] = S_mix[max_idx]
    end

    return nothing

end
"Calculate the Kolmogorov length scale `η` from viscousity and average TKE dissapation."
η(ν, ϵ) = (ν^3 / ϵ)^(1/4)
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

        NCDataset(file, "a") do ds
            Sc = ds.attrib["Sc"]
            find_num = findfirst(' ', ds.attrib["ν"]) - 1
            ν = parse(Float64, ds.attrib["ν"][1:find_num])
            η_min = (ν^3 / maximum(ds["ϵ_maximum"][:]))^(1/4)
            ds.attrib["η (min)"] = η_min # minimum space and time Kolmogorov scale
            ds.attrib["λ_B"] = η_min / sqrt(Sc) # minimum space and time Batchelor scale
        end

    elseif isequal(file_type, ".jld2")

        ϵ_ts = FieldTimeSeries(file, "ϵ_maximum", backend = OnDisk())
        Sc = load(file, "Non_dimensional_numbers")["Sc"]
        ν = load(file, "closure/ν")
        η_min = (ν^3 / maximum(ϵ_ts))^(1/4)

        jldopen(file, "a+") do f
            f["minimum_kolmogorov_scale"] = η_min
            f["minimum_batchelor_scale"] = η_min / sqrt(Sc)
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

    if simulation.output_writers[:tracers] isa NetCDFOutputWriter

        for key ∈ simulation.output_writers.keys
            if key != :checkpointer
                NCDataset(simulation.output_writers[key].filepath, "a") do ds
                    ds.attrib["EOS"] = summary(model.buoyancy.model.equation_of_state.seawater_polynomial)
                    ds.attrib["Reference density"] = "$(model.buoyancy.model.equation_of_state.reference_density)kgm⁻³"
                    ds.attrib["ν"]  = "$(model.closure.ν) m²s⁻¹"
                    ds.attrib["κₛ"] = "$(model.closure.κ.S) m²s⁻¹"
                    ds.attrib["κₜ"] = "$(model.closure.κ.T) m²s⁻¹"
                    for key ∈ keys(nd_nums)
                        ds.attrib[key] = nd_nums[key]
                    end
                end
            end
        end

    else

        #TODO: save this to all output writers
        jldopen(simulation.output_writers[:computed_output].filepath, "a+") do file
            file["Non_dimensional_numbers"] = nd_nums
        end

    end

    return nothing

end
"""
    function inferred_vertical_diffusivity(saved_output, flux, gradient)
Calculate the inferred vertical diffusivity using the horizontally averaged vertical
buoyancy gradient and horizontally averaged vertical buoyancy flux
```math
κᵥ = -\\frac{\\overline{w'b'}}{\\overline{\\frac{∂b}{∂z}}}.
```
**Note:** the horizontally avergaed vertical buoyancy gradient and vertical buoyancy flux
must be saved in `saved_output`. This is the default behaviour as of version 0.4.5.
Further it is assumed the horizontal resolution is equal.
"""
function inferred_vertical_diffusivity!(saved_output::AbstractString, flux::Symbol, vertical_gradient::Symbol)

    NCDataset(saved_output, "a") do ds
        vertical_gradient = ds[vertical_gradient][2:end, :]
        replace!(vertical_gradient, 0 => NaN)
        flux = ds[flux][:, :]
        ∫ₐκᵥ = similar(flux)
        ∫ₐκᵥ .= flux ./ vertical_gradient
        defVar(ds, "∫ₐκᵥ", ∫ₐκᵥ, ("zC", "time"),
               attrib = Dict("longname" => "Horizontally integrated inferred vertical diffusivity",
                             "units" => "m²s⁻¹"))
        dV = (diff(ds[:xC][1:2]) .* diff(ds[:yC][1:2])) .* diff(ds[:zF][:])
        replace!(∫ₐκᵥ, NaN => 0)
        ∫κᵥ = mapslices(sum, ∫ₐκᵥ .* dV, dims = 1)
        defVar(ds, "∫κᵥ", ∫κᵥ, ("time",),
        attrib = Dict("longname" => "Volume integrated inferred vertical diffusivity",
                      "units" => "m²s⁻¹"))
    end

    return nothing

end
"""
    function horizontal_average_profile!(tracers::AbstractString)
Compute and save the horizontally averaged salinity and temperature profiles that are
saved in `tracers`.
"""
function horizontal_average_profile!(tracers::AbstractString)

    NCDataset(tracers, "a") do ds

        z = ds[:zC]
        time = ds[:time]
        S, T = ds[:S], ds[:T]
        S_profile = Array{Float64}(undef, length(z), length(time))
        T_profile = Array{Float64}(undef, length(z), length(time))

        for t ∈ eachindex(time)

            S_profile[:, t] = reshape(mean(S[:, :, :, t], dims = (1, 2)), :)
            T_profile[:, t] = reshape(mean(T[:, :, :, t], dims = (1, 2)), :)

        end

        defVar(ds, "S_ha_profile", S_profile, ("zC", "time"),
               attrib = Dict("longname" => "Horizontally averaged salinity profile",
                             "units" => "gkg⁻¹"))
        defVar(ds, "T_ha_profile", T_profile, ("zC", "time"),
               attrib = Dict("longname" => "Horizontally averaged temperature profile",
                             "units" => "°C"))

    end

    return nothing
end
"""
    funciton find_file_type(file::AbstractString)
Return the file type (either `.nc` or `.jld2`) of a `file`.
"""
find_file_type(file::AbstractString) = file[findlast('.', file):end]
