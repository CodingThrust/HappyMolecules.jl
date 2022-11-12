using HappyMolecules
using Test, Random

using StaticArrays
using PythonCall
plt = pyimport("matplotlib.pyplot")

@testset "location initialization" begin
    Random.seed!(2)
    # random locations in an integer box size
    box = PeriodicBox(10, 10)
    @test HappyMolecules.largest_distance(box) == 5 * sqrt(2)
    locs = random_locations(box, 10000)
    @test length(locs) == 10000
    @test isapprox(sum(locs)/10000, SVector(4.5, 4.5); atol=5e-2)
    @test HappyMolecules.volume(box) == 100
    
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

@testset "binning" begin
    # constructor
    bin = Bin(-1.0, 1.0, 20)
    # ticks
    @test ticks(bin) ≈ [-1.05 + 0.1 * i for i=1:20]

    # push!
    @test_throws AssertionError push!(bin, -2.0)
    @test_throws AssertionError push!(bin, 2.0)
    r1 = zeros(Int, 20); r1[end] += 1
    @test push!(bin, 0.99).counts == r1
    r1[end] += 1
    @test push!(bin, 0.91).counts == r1
    r1[end-1] += 1
    @test push!(bin, 0.89).counts == r1

    # ncounts
    @test ncounts(bin) == 3

    # empty!
    empty!(bin)
    @test ncounts(bin) == 0
end

# Case study 4, close to the triple point of a Lennard-Jones Fluid
natoms = 108
temperature = 0.728
density = 0.8442
Nt = 2000
Δt = 0.001

# the box
volume = natoms / density
L = volume ^ (1/3)
box = PeriodicBox(SVector(L, L, L))

# initial status
lattice_pos = uniform_locations(box, natoms)
velocities = [rand(SVector{3, Float64}) .- 0.5 for _ = 1:natoms]
rc = L/2
md = molecule_dynamics(; lattice_pos, velocities, box, temperature, rc, Δt, compute_fr=true)

# Q: how to match the initial potential energy?
# Anderson thermalstat.
# Nose-Hoover thermalstat, difficult but better.
ps = Float64[]
ks = Float64[]
temps = Float64[]
for j=1:Nt
    step!(md)
    push!(ps, potential_energy(md))
    push!(ks, kinetic_energy(md))
    push!(temps, HappyMolecules.temperature(md))
end

# energy conservation
@test isapprox(ps[1] + ks[1], ps[end] + ks[end]; atol=1e-2)

fig = plt.figure(; figsize=(8, 6))
plt.plot(1:Nt, ps; label="Potential energy")
plt.plot(1:Nt, ks; label="Kinetic energy")
plt.plot(1:Nt, ps .+ ks; label="Total energy", color="k", ls="--")
plt.xlabel("step")
plt.ylabel("Energy/N")
plt.legend()
plt.show()

# ╔═╡ 0b7b46eb-8b33-48fa-b4e7-60003003cd3a
let
    filename = tempname() * ".mp4"
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
