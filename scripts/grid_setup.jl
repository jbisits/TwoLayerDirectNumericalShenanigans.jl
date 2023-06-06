## Setting up a variable spaced grid, with very high resolution over the middle 0.5m of the
#  domian and increasing/decreasing resolution at the upper/lower part of domain.
#  This is NOT WORKING YET.

resolution = (Nx = 10, Ny = 10 , Nz = 5000)

h(k) = -(k - 1) / resolution.Nz
function z_face_spacing(k)

    if h(k) ≥ -0.25
        tanh((h(k) + 0.26)/resolution.Nz)
    elseif h(k) ≤ -0.75
        -tanh((h(k) + 0.74)/resolution.Nz)
    else
        1e-5
    end

end
Δz_mm, Δz_cm = 1e-4, 1e-2

cumsum(z_face_spacing.(1:resolution.Nz + 1))
gaussian_spacing(k) = exp(-(((k - resolution.Nz/2)/ resolution.Nz )- 0.5)^2/2)
grid = RectilinearGrid(topology = (Periodic, Periodic, Bounded),
                        size = (resolution.Nx, resolution.Ny, resolution.Nz),
                        x = (-0.5/2, 0.5/2),
                        y = (-0.5/2, 0.5/2),
                        z = z_face_spacing.(1:resolution.Nz + 1))
scatterlines(zspacings(grid, Center()), znodes(grid, Center()))

## Two distinct spacings works fine
Δz_mm, Δz_cm = 1e-4, 1e-2
z_spacing_vector = vcat(-Lz:Δz_cm:-Lz /2 - 0.06,
                        -Lz /2 - 0.05:Δz_mm:-Lz /2 + 0.05,
                        -Lz /2 + 0.06:Δz_cm:0)

function cubically_spaced_faces(k)
    if  k < resolution.Nz / 2
        -((1 - k)/(1 - resolution.Nz/2))^(1/2)/2 - 0.5
    else k > resolution.Nz / 2
        ((resolution.Nz + 1 - k)/(1 + resolution.Nz/2))^(1/2)/2 - 0.5
    #else

    end
end
sinh_spaced_faces(k) = sinh((k - resolution.Nz/2) / (resolution.Nz + 1)) - 0.5
function tanh_spaced_faces(k)
    if k ≤ resolution.Nz / 4
        -1
        #1/2 - tanh((k - resolution.Nz) / (1 - resolution.Nz / 4)) / 2tanh(1)
    elseif k ≥ 3*resolution.Nz / 4
        0
        #1/2 + tanh((k - 3*resolution.Nz / 4) / (resolution.Nz + 1 - 3*resolution.Nz / 4 )) / 2tanh(1)
    else
        #tanh((k - resolution.Nz/4) / (resolution.Nz / 4))
        1
    end
end
σ = 2
hyperbolically_spaced_faces(k) = -1 * (1 - tanh(σ * (k - 1) / resolution.Nz)/ tanh(σ))
refinement = 1.8 # controls spacing near surface (higher means finer spaced)
stretching = 100  # controls rate of stretching at bottom

# Normalized height ranging from 0 to 1
h_(k) = (k - 1) / resolution.Nz

# Linear near-surface generator
ζ₀(k) = 1 + (h_(k) - 1) / refinement

# Bottom-intensified stretching function
Σ(k) = (1 - exp(-stretching * h_(k))) / (1 - exp(-stretching))

# Generating function
z_faces(k) = 1 * (ζ₀(k) * Σ(k) - 1)
grid = RectilinearGrid(topology = (Periodic, Periodic, Bounded),
                        size = (resolution.Nx, resolution.Ny, resolution.Nz),
                        x = (-0.5/2, 0.5/2),
                        y = (-0.5/2, 0.5/2),
                        z = z_faces)
fig, ax = scatterlines(zspacings(grid, Center()), znodes(grid, Center()))
ax.title = "Grid spacing for DNS"
ax.xlabel = "Δz"
ax.ylabel = "z"
fig
