module HappyMolecules

using StaticArrays
using Statistics
using DocStringExtensions

export molecule_dynamics, step!
export PeriodicBox, Box, random_locations, uniform_locations, volume
export positions, velocities, forces, num_particles,
        kinetic_energy, temperature, potential_energy,
        pressure

include("Core.jl")

end
