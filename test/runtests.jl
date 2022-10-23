using HappyMolecules
using Test

using StaticArrays
natoms = 100
box = PeriodicBox(SVector(20.0, 20.0, 20.0))
x0 = random_locations(box, natoms)
rc2 = 10.0
temperature = 1.0
md = molecule_dynamics(x0, box, temperature, rc2, 0.01)
step!(md);

@testset "HappyMolecules.jl" begin
    # Write your tests here.
end
