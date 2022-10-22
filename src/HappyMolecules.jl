module HappyMolecures

function md(lattice_pos, temperature, rc2)
    ############# INIT ###############
    npart = length(lattice_pos)
    # initialize locations as x and velocities as v
    x = copy(lattice_pos)
    v = rand(npart) .- 0.5

    # since we have degree of freedoms 3
    # m*v^2/2 = 3/2*k*T
    # Q: why the rescaling factor is not computed after subtracting the velocity center?
    fs = sqrt(3 * temperature / sum(abs2, v))

    v_mean = sum(v) / npart
    # set velocity center of mass to zero.
    v .= (v .- v_mean) .* fs
    # position previous time step
    xm = x .- v .* t

    # assert rc < box / 2
    t = 0
    while t < tmax
        # compute the force
        field, energy = force!(field, x, v, rc2, box)
        integrate!(x, xm, field, energy)
        t = t + δt
    end
end

# Lennard-Jones potential
# f(r) = 48 r⃗ / r² (1 / r¹² - 0.5* 1 / r⁶)
function force!(field::AbstractVector{T}, x::AbstractVector, v::AbstractVector, rc2, box::Box) where T
    npart = length(x)
    @debug @assert length(v) == length(field) == npart
    fill!(field, zero(T))
    energy = zero(T)

    for i=1:npart, j=i+1:npart
        xr = distance_vector(x[i], x[j], box)
        xr2 = sum(abs2, xr)
        if xr2 < rc2
            r2i = 1 / r2
            r6i = r2 ^ 3

            # Lennard-Jones potential
            ff = 48 * r2i * r6i * (r6i - 0.5)
            field[i] .+= ff .* xr
            field[j] .-= ff .* xr

            # update energy
            energy += 4 * r6i * (r6i - 1) - ecut
        end
    end
    return field, energy
end

function distance_vector(x, y, box::PeriodicBox)
    r = x - y
    return round.(r ./ box.dimensions) .* box.dimensions
end

function integrate!(x, xm, field, energy)
    sumv = 0
    sumv2 = 0

    # the Verlet algorithm
    for i=1:npart
        xx  = 2 * x[i] - xm[i] + δ^2 * field[i]
        vi = (xx - xm[i]) / 2δ
        sumv += vi
        sumv2 = sumv2 + sum(abs2, vi)
        xm[i] = x[i]
        x[i] = xx
    end
    temperature = sumv2 / (3 * npart)
    energy_total = (energy + 0.5 * sumv2) / npart
    return temperature, energy_total
end

end
