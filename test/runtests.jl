using HappyMolecules
using Test

using StaticArrays
natoms = 100
x0 = rand(SVector{3, Float64}, natoms)
box = PeriodicBox(SVector(1.0, 1.0, 1.0))
rc2 = 0.1
temperature = 1.0
md = molecule_dynamics(x0, box, temperature, rc2, 0.01)
step!(md)

@testset "HappyMolecules.jl" begin
    # Write your tests here.
end
