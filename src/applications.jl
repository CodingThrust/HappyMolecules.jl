module Applications

using Random, StaticArrays
using ..HappyMolecules

export lennard_jones_triple_point

# Case study 4, close to the triple point of a Lennard-Jones Fluid in a 3D periodic box
function lennard_jones_triple_point(;
        natoms::Int = 108,
        temperature::Real = 0.728,
        density::Real = 0.8442,
        Nt = 2000,
        Δt = 0.001,
        seed::Int = 2,
        gr_lastn::Int = 500,
    )
    Random.seed!(seed)

    # the box
    volume = natoms / density
    L = volume ^ (1/3)
    box = PeriodicBox(SVector(L, L, L))

    # initial status
    lattice_pos = uniform_locations(box, natoms)
    velocities = [rand(SVector{3, Float64}) .- 0.5 for _ = 1:natoms]
    rc = L/2

    # create a `MDRuntime` instance
    md = molecule_dynamics(; lattice_pos, velocities, box, temperature, rc, Δt, potential=LennardJones(; rc))

    # Q: how to match the initial potential energy?
    # Anderson thermalstat.
    # Nose-Hoover thermalstat, difficult but better.
    ps = Float64[]
    ks = Float64[]
    temps = Float64[]

    bin = Bin(0.0, L/2, 200)
    for j=1:Nt
        step!(md)
        push!(ps, mean_potential_energy(md))
        push!(ks, mean_kinetic_energy(md))
        push!(temps, HappyMolecules.temperature(md))
        if j > Nt - gr_lastn
            HappyMolecules.collect_gr!(md, bin)
        end
    end
    return (;
        runtime = md,
        potential_energy = ps,
        kinetic_energy = ks,
        radial_ticks = ticks(bin),
        radial_distribution = HappyMolecules.finalize_gr(md, bin, gr_lastn)
    )
end

end
