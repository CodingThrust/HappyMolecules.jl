using HappyMolecules
using Test, Random

using StaticArrays

@testset "location initialization" begin
    Random.seed!(2)
    # random locations in an integer box size
    box = PeriodicBox(10, 10)
    locs = random_locations(box, 10000)
    @test length(locs) == 10000
    @test isapprox(sum(locs)/10000, SVector(4.5, 4.5); atol=5e-2)
    
    # random locations in a floating point box size
    box = PeriodicBox(10.0, 10.0)
    locs = random_locations(box, 10000)
    @test length(locs) == 10000
    @test isapprox(sum(locs)/10000, SVector(5.0, 5.0); atol=5e-2)

    # uniform locations in a floating point box size
    box = PeriodicBox(10.0, 10.0)
    locs = uniform_locations(box, 4)
    @test length(locs) == 4
    @test locs ≈ [SVector(0.0, 0.0), SVector(5.0, 0.0), SVector(0.0, 5.0), SVector(5.0, 5.0)]
end

# Case study 4, close to the triple point of a Lennard-Jones Fluid
natoms = 108
temperature = 0.728
density = 0.8442
filename = tempname() * ".mp4"

# the box
volume = natoms / density
L = volume ^ (1/3)
box = PeriodicBox(SVector(L, L, L))

# initial status
lattice_pos = uniform_locations(box, natoms)
velocities = [rand(SVector{3, Float64}) .- 0.5 for _ = 1:natoms]
rc = L/2
md = molecule_dynamics(; lattice_pos, velocities, box, temperature, rc, Δt=0.0005)
step!(md);

using GLMakie

# ╔═╡ 0b7b46eb-8b33-48fa-b4e7-60003003cd3a
let
    fig = Figure(; resolution=(800, 800))
    ax = Axis3(fig[1,1]; aspect=:data)
    limits = GLMakie.FRect3D((0, 0, 0),(box.dimensions...,))
    limits!(ax, limits)
	points = Observable([Point3f(x...,) for x in md.x])
	directions = Observable([Point3f(x/100...,) for x in md.field])
	scatter!(ax, points)
	arrows!(ax, points, directions; linewidth=0.02, arrowsize=0.1)
	record(fig, filename, 1:500; framerate = 30, sleep=true) do i
		for j=1:10
			step!(md)
		end
		points[] = [Point3f(mod.(x, box.dimensions)...,) for x in md.x]
		directions[] = [Point3f(x/100...,) for x in md.field]
	end
end

@testset "HappyMolecules.jl" begin
    # Write your tests here.
end
