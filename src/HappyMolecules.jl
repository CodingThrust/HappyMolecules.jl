module HappyMolecules

using StaticArrays
using Statistics

export molecule_dynamics, step!
export PeriodicBox, Box, random_locations, uniform_locations

abstract type Box{D} end

struct PeriodicBox{D, T} <: Box{D}
    dimensions::SVector{D,T}
end
function PeriodicBox(x::T, xs::T...) where T<:Real
    return PeriodicBox(SVector(x, xs...))
end

function random_locations(box::PeriodicBox{D, T}, natoms::Int) where {D,T}
    return [rand(SVector{D, T}) .* box.dimensions for _=1:natoms]
end

function uniform_locations(box::PeriodicBox{D, T}, natoms::Int) where {D,T}
    L = round(Int, natoms ^ (1/D))
    L^D ≈ natoms || error("`natoms` should be L^$D for some integer L.")
    return vec([SVector((idx.I .- 1) ./ L .* box.dimensions) for idx in CartesianIndices(ntuple(i->L, D))])
end

norm2(v::SVector) = sum(abs2, v)

Base.@kwdef struct MDConfig{D, BT<:Box{D}, RT}
    box::BT
    n::Int
    temperature::RT
    rc2::RT
    Δt::RT
end

mutable struct MDRuntime{D, BT, T}
    const config::MDConfig{D,BT,T}
    t::T
    energy::T
    const x::AbstractVector{SVector{D, T}}
    const xm::AbstractVector{SVector{D, T}}
    const field::AbstractVector{SVector{D, T}}
end

function molecule_dynamics(lattice_pos::AbstractVector{SVector{D, T}}, box::Box{D}, temperature::Real, rc2::Real, Δt::Real) where {D, T}
    # assert rc < box / 2
    ############# INIT ###############
    n = length(lattice_pos)
    # initialize locations as x and velocities as v
    x = copy(lattice_pos)
    v = [rand(SVector{D, T}) .- 0.5 for _ = 1:n]

    # since we have degree of freedoms 3
    # m*v^2/2 = D/2*k*T, because we have `D` degrees of freedoms to move.
    # Q: why the rescaling factor is not computed after subtracting the velocity center?
    fs = sqrt(D * temperature / mean(norm2, v))

    v_mean = mean(v)
    # set velocity center of mass to zero.
    v .= (v .- Ref(v_mean)) .* fs
    @debug "mean v² = $(mean(norm2, v))"
    # position previous time step
    xm = x .- v .* Δt

    # intialize a vector field
    config = MDConfig(; box, n, temperature, rc2, Δt)
    return MDRuntime(config, zero(T), zero(T), x, xm, zero(v))
end

function step!(md::MDRuntime)
    # compute the force
    md.energy = force!(md.field, md.x, md.config.rc2, md.config.box)
    integrate!(md.x, md.xm, md.field, md.energy, md.config.Δt)
    md.t += md.config.Δt
    return md
end

# Lennard-Jones potential
# f(r) = 48 r⃗ / r² (1 / r¹² - 0.5* 1 / r⁶)
function force!(field::AbstractVector{SVector{D, T}}, x::AbstractVector{SVector{D, T}}, rc2, box::Box) where {D,T}
    npart = length(x)
    @debug @assert length(v) == length(field) == npart
    fill!(field, zero(SVector{D, T}))
    energy = zero(T)
    ecut = 4 * (1/rc2^6 - 1/rc2^3)

    for i=1:npart, j=i+1:npart
        xr = distance_vector(x[i], x[j], box)
        r2 = norm2(xr)
        if r2 < rc2
            r2i = 1 / r2
            r6i = r2i ^ 3

            # Lennard-Jones potential
            ff = 48 * r2i * r6i * (r6i - 0.5)
            field[i] += ff * xr
            field[j] -= ff * xr

            # update energy
            # Q: why minus the ecut?
            energy += 4 * r6i * (r6i - 1) - ecut
        end
    end
    @debug "mean energy = $(energy / npart), ecut = $ecut"
    return energy
end

function distance_vector(x, y, box::PeriodicBox)
    r = x - y
    return r .- round.(r ./ box.dimensions) .* box.dimensions
end

function integrate!(x::AbstractVector{SVector{D,T}}, xm, field, energy, Δt) where {D,T}
    sumv = zero(SVector{D,T})
    sumv2 = zero(T)
    npart = length(x)

    # the Verlet algorithm
    for i=1:npart
        xx  = 2 * x[i] - xm[i] + Δt^2 * field[i]
        vi = (xx - xm[i]) / (2 * Δt)
        sumv += vi
        sumv2 += norm2(vi)
        xm[i] = x[i]
        x[i] = xx
    end
    temperature = sumv2 / (D * npart)
    energy_total = (energy + 0.5 * sumv2) / npart
    return temperature, energy_total
end

end
