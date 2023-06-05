## Setting up a variable spaced grid, with very high resolution over the middle 0.5m of the
#  domian and increasing/decreasing resolution at the upper/lower part of domain.
#  This is NOT WORKING YET.

resolution = (Nx = 10, Ny = 10 , Nz = 100)

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
gaussian_spacing(k) = exp(-((k - 1)/ resolution.Nz - 0.5)^2/2)
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

cubically_spaced_faces(k) = k < resolution.Nz / 2 ?
                            -((k - resolution.Nz/2)/(1 - resolution.Nz/2))^(1/2)/2 - 0.5 :
                            ((k - resolution.Nz/2)/(resolution.Nz/2))^(1/2)/2 - 0.5
sinh_spaced_faces(k) = sinh((k - resolution.Nz/2) / (resolution.Nz + 1)) - 0.5
tanh_spaced_faces(k) = k < resolution.Nz / 2 ?
                        -cosh((k - resolution.Nz/2) / (1 + resolution.Nz)) -0.5 :
                        cosh((k - resolution.Nz/2) / (1 + resolution.Nz)) -0.5
grid = RectilinearGrid(topology = (Periodic, Periodic, Bounded),
                        size = (resolution.Nx, resolution.Ny, resolution.Nz),
                        x = (-0.5/2, 0.5/2),
                        y = (-0.5/2, 0.5/2),
                        z = tanh_spaced_faces)
scatterlines(zspacings(grid, Center()), znodes(grid, Center()))
