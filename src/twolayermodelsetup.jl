
"""
    function set_two_layer_initial_conditions(model::AbstractModel, S::NamedTuple, Θ::NamedTuple)
Set initial conditions for temperature and salinity in a two layer model. Initial values for
the absolute salinity (`S`) and conservative temperature (`Θ`) in each layer must be passed
to the function as a `NamedTuple` in the form `S = (upper = value, lower = value)` as well
as the model where the initial conditions should be set. **Note** the interface of the
layers is in the middle of the domain.
"""
function set_two_layer_initial_conditions!(model::Oceananigans.AbstractModel, S::NamedTuple,
                                          Θ::NamedTuple)

    S₀ = Array{Float64}(undef, size(model.grid))
    Θ₀ = Array{Float64}(undef, size(model.grid))

    for i ∈ axes(S₀, 3)

        if i ≤ floor(model.grid.Nz / 2)
            S₀[:, :, i] .= S.lower
            Θ₀[:, :, i] .= Θ.lower
        else
            S₀[:, :, i] .= S.upper
            Θ₀[:, :, i] .= Θ.upper
        end

    end

    set!(model, T = Θ₀, S = S₀)

    return nothing

end

"""
    function viscosity_and_diffusivity(Prandtl_number::Number, Raleigh_number::Number)
Return the kinematic viscosity (`ν`) and the diffusivity (`κ`) from the
[Prandtl number](https://en.wikipedia.org/wiki/Prandtl_number) and
the [Raleigh number](https://en.wikipedia.org/wiki/Rayleigh_number).
The Prandtl number is
```math
    Pr = \\frac{ν}{κ}
```
and the Raleigh number is
```math
    Ra = \\frac{ΔρgL³}{μκ} = \\frac{bL³}{νκ}.
```
Using these expressions we can derive `ν` and `κ` in terms of the non-dimensional numbers,
the buoyancy and the domain size
```math
ν = \\sqrt{\\frac{PrbL³}{Ra}} \\quad κ = \\sqrt{\\frac{bL³}{RaPr}}.
```
"""
function viscosity_and_diffusivity(Prandtl_number::Number, Raleigh_number::Number)

    return nothing
end
