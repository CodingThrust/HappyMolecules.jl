module HappyMolecules

using StaticArrays
using Statistics
using DocStringExtensions

export molecule_dynamics, step!
export PeriodicBox, Box, random_locations, uniform_locations
export positions, velocities, forces, num_particles, kinetic_energy, temperature

include("Core.jl")

end
