abstract type Box{D} end

struct PeriodicBox{D, T} <: Box{D}
    dimensions::SVector{D,T}
end
function PeriodicBox(x::T, xs::T...) where T<:Real
    return PeriodicBox(SVector(x, xs...))
end

function random_locations(box::PeriodicBox{D, T}, natoms::Int) where {D,T}
    return [rand_uniform.(Ref(zero(T)), box.dimensions) for _=1:natoms]
end
function rand_uniform(min::T, max::T) where T <: AbstractFloat
    return rand(T) * (max - min) + min
end
function rand_uniform(min::T, max::T) where T <: Integer
    return rand(0:max-min-1) + min
end

function uniform_locations(box::PeriodicBox{D, T}, natoms::Int) where {D,T}
    L = ceil(Int, natoms ^ (1/D))
    L^D ≈ natoms || @warn("`natoms = $natoms` is not equal to L^$D for some integer L.")
    CIS = CartesianIndices(ntuple(i->L, D))
    return vec([SVector((CIS[i].I .- 1) ./ L .* box.dimensions) for i=1:natoms])
end

norm2(v::SVector) = sum(abs2, v)

Base.@kwdef struct MDConfig{D, BT<:Box{D}, RT}
    box::BT
    n::Int
    temperature::RT
    rc2::RT
    Δt::RT
    compute_fr::Bool
end

mutable struct MDRuntime{D, BT, T}
    const config::MDConfig{D,BT,T}
    t::T
    potential_energy::T
    mean_fr::T
    const x::Vector{SVector{D, T}}
    const xm::Vector{SVector{D, T}}
    const v::Vector{SVector{D,T}}
    const field::Vector{SVector{D, T}}
end

# get properties from the run time information.
positions(r::MDRuntime) = r.x
velocities(r::MDRuntime) = r.v
forces(r::MDRuntime) = r.field
num_particles(r::MDRuntime) = r.config.n
kinetic_energy(r::MDRuntime{D}) where D = temperature(r) * (D / 2)
potential_energy(r::MDRuntime{D}) where D = r.potential_energy / r.config.n
mean_fr(r::MDRuntime{D}) where D = r.mean_fr

"""
$TYPEDSIGNATURES

The temperature ``T`` is measured by computing the average kinetic energy per degree of freedom.
```math
k_B T = \\frac{\\langle 2 \\mathcal{K} \\rangle}{f}.
```
"""
function temperature(r::MDRuntime{D, BT, T}) where {D, BT, T}
    npart = num_particles(r)
    sumv2 = zero(T)
    for i = 1:npart
        vi = r.v[i]
        sumv2 += norm2(vi)
    end
    return sumv2 / (D * npart)
end

"""
Compute the constant volume capacity is using the following equation
```math
\\langle K^2 \\rangle - \\langle K \\rangle^2 = \\frac{3 k_b^2 T^2}{2N}(1-\\frac{3k_B}{2C_v})
```
"""
function heat_capacity(r::MDRuntime{D}) where D
    npart = num_particles(r)
    sumv2 = zero(T)
    sumv = zero(SVector{D,T})
    for i = 1:npart
        vi = r.v[i]
        sumv2 += norm2(vi)
        sumv += vi
    end
    fluctuation = sumv2 / npart - (sumv/npart) ^ 2
    t2 = 3 * temperature ^ 2 / 2 / num_particles(r)
    t3 = (1 - t2 / fluctuation)  # = 3/(2Cv)
    return 1.5 / t3
end

function radial_distribution()
end

"""
The most common among the ways to measure the pressure ``P`` is based on the virial equation for the pressure.
```math
P = \\rho k_B T + \\frac{1}{dV}\\langle\\sum_{i<j} f(r_{ij}) \\cdot r_{ij}\\rangle
```
"""
pressure(r::MDRuntime{D}) where D = pressure_formula(density(r), temperature(r), mean_fr(r), D, volume(r))
function pressure_formula(ρ, T, mean_fr, D, volume)
    ρ * T + mean_fr / D/ volume
end

function molecule_dynamics(; lattice_pos::AbstractVector{SVector{D, T}}, velocities::AbstractVector{SVector{D, T}}, box::Box{D}, temperature::Real, rc::Real, Δt::Real, compute_fr::Bool=false) where {D, T}
    # assert rc < box / 2
    ############# INIT ###############
    n = length(lattice_pos)
    # initialize locations as x and velocities as v
    x = copy(lattice_pos)
    v = copy(velocities)

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
    config = MDConfig(; box, n, temperature, rc2=rc^2, Δt, compute_fr)
    return MDRuntime(config, zero(T), zero(T), zero(T), x, xm, v, zero(v))
end

function step!(md::MDRuntime)
    # compute the force
    md.potential_energy, md.mean_fr = force!(md.field, md.x, md.config.rc2, md.config.box; compute_fr=md.config.compute_fr)
    integrate!(md.x, md.xm, md.v, md.field, md.config.Δt)
    md.t += md.config.Δt
    return md
end

# Lennard-Jones potential
# f(r) = 48 r⃗ / r² (1 / r¹² - 0.5* 1 / r⁶)
function force!(field::AbstractVector{SVector{D, T}}, x::AbstractVector{SVector{D, T}}, rc2, box::Box; compute_fr::Bool) where {D,T}
    npart = length(x)
    @debug @assert length(v) == length(field) == npart
    fill!(field, zero(SVector{D, T}))
    potential_energy = zero(T)
    ecut = 4 * (1/rc2^6 - 1/rc2^3)
    fr = zero(T)

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

            # compute ⟨f⃗ ⋅ r⃗⟩ e.g. for computing the pressure
            if compute_fr
                fr += ff * xr^2
            end

            # update potential_energy
            # Q: why minus the ecut?
            potential_energy += 4 * r6i * (r6i - 1) - ecut
        end
    end
    @debug "mean potential energy = $(potential_energy / npart), ecut = $ecut"
    return potential_energy, fr / (npart * (npart-1) ÷ 2)
end

function distance_vector(x, y, box::PeriodicBox)
    r = x - y
    return r .- round.(r ./ box.dimensions) .* box.dimensions
end

function integrate!(x::AbstractVector{SVector{D,T}}, xm, v, field, Δt) where {D,T}
    npart = length(x)

    # the Verlet algorithm
    for i=1:npart
        xx  = 2 * x[i] - xm[i] + Δt^2 * field[i]
        v[i] = (xx - xm[i]) / (2 * Δt)
        xm[i] = x[i]
        x[i] = xx
    end
end

