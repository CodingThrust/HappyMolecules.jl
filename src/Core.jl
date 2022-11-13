abstract type PotentialField end
abstract type Box{D} end

struct PeriodicBox{D, T} <: Box{D}
    dimensions::SVector{D,T}
end
function PeriodicBox(x::T, xs::T...) where T<:Real
    return PeriodicBox(SVector(x, xs...))
end
volume(box::PeriodicBox) = prod(box.dimensions)

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
largest_distance(box::PeriodicBox) = sqrt(sum(abs2, box.dimensions)) / 2

norm2(v::SVector) = sum(abs2, v)

Base.@kwdef struct MDConfig{D, RT, BT<:Box{D}, PT<:PotentialField}
    box::BT
    potential::PT
    n::Int
    temperature::RT
    rc2::RT
    Δt::RT
end

mutable struct MDRuntime{D, T, BT, PT}
    const config::MDConfig{D,T,BT,PT}
    t::T
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
function sum_fr(md::MDRuntime{D}) where D
    # compute ⟨f⃗ ⋅ r⃗⟩ e.g. for computing the pressure
    npart = num_particles(md)
    fr = zero(T)
    for i=1:npart-1
        for j=i+1:npart
            xr = distance_vector(md.x[i], md.x[j], md.config.box)
            fr += sum(force(md.config.potential, xr) .* xr)
        end
    end
    return fr
end

mean_kinetic_energy(r::MDRuntime{D}) where D = temperature(r) * (D / 2)
"""
$TYPEDSIGNATURES
"""
function mean_potential_energy(r::MDRuntime{D,T}) where {D,T}
    npart = num_particles(r)
    eng = zero(T)
    for i=1:npart-1
        for j=i+1:npart
            eng += potential_energy(r.config.potential, sqrt(norm2(r.x[i] - r.x[j])))
        end
    end
    @debug "mean potential energy = $(eng / npart)"
    return eng / npart
end

"""
$TYPEDSIGNATURES

The temperature ``T`` is measured by computing the average kinetic energy per degree of freedom.
```math
k_B T = \\frac{\\langle 2 \\mathcal{K} \\rangle}{f}.
```
"""
function temperature(r::MDRuntime{D, T}) where {D, T}
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

"""
The most common among the ways to measure the pressure ``P`` is based on the virial equation for the pressure.
```math
P = \\rho k_B T + \\frac{1}{dV}\\langle\\sum_{i<j} f(r_{ij}) \\cdot r_{ij}\\rangle
```
"""
pressure(r::MDRuntime{D}) where D = pressure_formula(density(r), temperature(r), sum_fr(r), D, volume(r.config.box))
function pressure_formula(ρ, T, sum_fr, D, volume)
    ρ * T + sum_fr / D/ volume
end
density(md::MDRuntime{D}) where D = md.config.n / volume(md.config.box)

function molecule_dynamics(; lattice_pos::AbstractVector{SVector{D, T}}, velocities::AbstractVector{SVector{D, T}}, box::Box{D}, potential::PotentialField, temperature::Real, rc::Real, Δt::Real) where {D, T}
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
    config = MDConfig(; box, n, temperature, rc2=rc^2, Δt, potential)
    return MDRuntime(config, zero(T), x, xm, v, zero(v))
end

function step!(md::MDRuntime)
    # compute the force
    update_force_field!(md.config.potential, md.field, md.x, md.config.rc2, md.config.box)
    integrate!(md.x, md.xm, md.v, md.field, md.config.Δt)
    md.t += md.config.Δt
    return md
end

# Benhcmark: lammps, gromacs
struct Bin{T}
    counts::Vector{Int}
    min::T
    max::T
end
function Bin(min::T, max::T, n::Int) where T
    counts = zeros(Int, n)
    return Bin(counts, min, max)
end

function ticks(bin::Bin)
    nticks = length(bin.counts)
    step = (bin.max - bin.min) / nticks
    return [bin.min + step * (i - 0.5) for i=1:nticks]
end
function Base.push!(bin::Bin{T}, val::T) where T
    @assert val >= bin.min && val < bin.max "the value $val is out of binning range: [$(bin.min), $(bin.max))."
    nticks = length(bin.counts)
    step = (bin.max - bin.min) / nticks
    bin.counts[floor(Int, (val - bin.min) / step) + 1] += 1
    return bin
end
function Base.empty!(bin::Bin{T}) where T
    bin.counts .= 0
    return bin
end
ncounts(bin::Bin) = sum(bin.counts)

function measure_gr(md::MDRuntime{D}; nbins=500, min_distance=0.0, max_distance=minimum(md.config.box.dimensions)/2) where D
    bin = Bin(min_distance, max_distance, nbins)
    measure_gr!(md, bin)
    finalize_gr(md, bin, 1)
end

# normalize over volume
function finalize_gr(md::MDRuntime, bin::Bin, niters::Int)
    nbins = length(bin.counts)
    Δr = (bin.max - bin.min) / nbins
    rs = ticks(bin)
    return map(1:nbins) do i
        num_particles = (4/3*π*((rs[i] + Δr/2)^3 - (rs[i] - Δr/2)^3) * density(md))
        (bin.counts[i] * 2) / niters / md.config.n / num_particles
    end
end

function collect_gr!(md::MDRuntime, bin::Bin)
    npart = md.config.n
    for i=1:npart-1, j=i+1:npart
        xr = distance_vector(md.x[i], md.x[j], md.config.box)
        r = sqrt(norm2(xr))
        if r < bin.max
            push!(bin, r)
        end
    end
    return bin
end

# Lennard-Jones potential
# f(r) = 48 r⃗ / r² (1 / r¹² - 0.5* 1 / r⁶)
function update_force_field!(potential::PotentialField, field::AbstractVector{SVector{D, T}}, x::AbstractVector{SVector{D, T}}, rc2, box::Box) where {D,T}
    npart = length(x)
    @assert length(field) == npart
    fill!(field, zero(SVector{D, T}))

    for i=1:npart, j=i+1:npart
        xr = distance_vector(x[i], x[j], box)
        r2 = norm2(xr)
        if r2 < rc2
            # Lennard-Jones potential
            f = force(potential, xr)
            field[i] += f
            field[j] -= f
        end
    end
end

"""
```math
f(r) = \\frac{48 \\vec{r}}{r^2}(\\frac{1}{r^{12}} - 0.5 \\frac{1}{r^6})
```
"""
struct LennardJones <: PotentialField
end

function force(::LennardJones, distance_vector::SVector)
    r2 = norm2(distance_vector)
    r2i = inv(r2)
    r6i = r2i ^ 3
    ff = 48 * r2i * r6i * (r6i - 0.5)
    return ff * distance_vector
end

function potential_energy(::LennardJones, distance::Real)
    r6i = distance ^ -6
    # Q: why minus ecut?
    rc2 = 2.519394287073761 ^ 2
     ecut = 4 * (1/rc2^6 - 1/rc2^3)
    return 4 * r6i * (r6i - 1) - ecut
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

