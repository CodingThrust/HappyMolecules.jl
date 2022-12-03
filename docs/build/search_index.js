var documenterSearchIndex = {"docs":
[{"location":"ref/","page":"Reference","title":"Reference","text":"Modules = [HappyMolecules, HappyMolecules.Applications]","category":"page"},{"location":"ref/#HappyMolecules.Bin","page":"Reference","title":"HappyMolecules.Bin","text":"The type for binning statistics.\n\n\n\n\n\n","category":"type"},{"location":"ref/#HappyMolecules.LennardJones","page":"Reference","title":"HappyMolecules.LennardJones","text":"f(r) = frac48 vecrr^2(frac1r^12 - 05 frac1r^6)\n\n\n\n\n\n","category":"type"},{"location":"ref/#HappyMolecules.heat_capacity-Union{Tuple{HappyMolecules.MDRuntime{D}}, Tuple{D}} where D","page":"Reference","title":"HappyMolecules.heat_capacity","text":"Compute the constant volume capacity is using the following equation\n\nlangle K^2 rangle - langle K rangle^2 = frac3 k_b^2 T^22N(1-frac3k_B2C_v)\n\nheat_capacity(md)\n\ndefined at /home/runner/work/HappyMolecules.jl/HappyMolecules.jl/src/Core.jl:119.\n\n\n\n\n\n","category":"method"},{"location":"ref/#HappyMolecules.measure_gr-Union{Tuple{HappyMolecules.MDRuntime{D}}, Tuple{D}} where D","page":"Reference","title":"HappyMolecules.measure_gr","text":"Measure the radial distribtion over a molecular dynamics runtime instance.\n\nArguments\n\nmd is the molecular dynamics runtime instance.\n\nKeyword argument\n\nnbins is the number of bins,\nmin_distance and max_distance are the minimum and maximum distance of the bins.\n\nmeasure_gr(md; nbins, min_distance, max_distance)\n\ndefined at /home/runner/work/HappyMolecules.jl/HappyMolecules.jl/src/Core.jl:245.\n\n\n\n\n\n","category":"method"},{"location":"ref/#HappyMolecules.pressure-Union{Tuple{HappyMolecules.MDRuntime{D}}, Tuple{D}} where D","page":"Reference","title":"HappyMolecules.pressure","text":"The most common among the ways to measure the pressure P is based on the virial equation for the pressure.\n\nP = rho k_B T + frac1dVlanglesum_ij f(r_ij) cdot r_ijrangle\n\npressure(md)\n\ndefined at /home/runner/work/HappyMolecules.jl/HappyMolecules.jl/src/Core.jl:140.\n\n\n\n\n\n","category":"method"},{"location":"ref/#HappyMolecules.random_locations-Union{Tuple{T}, Tuple{D}, Tuple{PeriodicBox{D, T}, Int64}} where {D, T}","page":"Reference","title":"HappyMolecules.random_locations","text":"random_locations(box::Box, natoms::Int) -> SVector\n\nReturns a set of random locations in a box\n\nrandom_locations(box, natoms)\n\ndefined at /home/runner/work/HappyMolecules.jl/HappyMolecules.jl/src/Core.jl:22.\n\n\n\n\n\n","category":"method"},{"location":"ref/#HappyMolecules.temperature-Union{Tuple{HappyMolecules.MDRuntime{D, T}}, Tuple{T}, Tuple{D}} where {D, T}","page":"Reference","title":"HappyMolecules.temperature","text":"The temperature T is measured by computing the average kinetic energy per degree of freedom.\n\nk_B T = fraclangle 2 mathcalK ranglef\n\ntemperature(md)\n\ndefined at /home/runner/work/HappyMolecules.jl/HappyMolecules.jl/src/Core.jl:103.\n\n\n\n\n\n","category":"method"},{"location":"ref/#HappyMolecules.ticks-Tuple{Bin}","page":"Reference","title":"HappyMolecules.ticks","text":"ticks(bin)\n\n\nReturn the ticks (center of boxes) of bins.\n\nticks(bin)\n\ndefined at /home/runner/work/HappyMolecules.jl/HappyMolecules.jl/src/Core.jl:211.\n\n\n\n\n\n","category":"method"},{"location":"ref/#HappyMolecules.uniform_locations-Union{Tuple{T}, Tuple{D}, Tuple{PeriodicBox{D, T}, Int64}} where {D, T}","page":"Reference","title":"HappyMolecules.uniform_locations","text":"uniform_locations(box, natoms)\n\n\nuniform_locations(box::Box, natoms::Int) -> SVector\n\nReturns a set of uniform locations in a box\n\nuniform_locations(box, natoms)\n\ndefined at /home/runner/work/HappyMolecules.jl/HappyMolecules.jl/src/Core.jl:39.\n\n\n\n\n\n","category":"method"},{"location":"ref/#HappyMolecules.Applications.lennard_jones_triple_point-Tuple{}","page":"Reference","title":"HappyMolecules.Applications.lennard_jones_triple_point","text":"lennard_jones_triple_point(\n;\n    natoms,\n    temperature,\n    density,\n    Nt,\n    Δt,\n    seed,\n    gr_lastn\n)\n\n\nCase study in Chapter 4 of the book \"Understanding Molecular Simulation, From Algorithms to Applications\". It is about the molecule dynamics simulation of a Lennard-Jones Fluid in a 3D periodic box. The parameters are set close to the triple point.\n\nKeyword arguments\n\nnatoms is the number of atoms.\ntemperature is the initial temperature.\ndensity is the density of atoms.\nNt is the number of tims steps.\nΔt is the time step.\nseed is the random seed.\ngr_lastn is the number of last n samples for collecting radial distribution.\n\nlennard_jones_triple_point(\n;\n    natoms,\n    temperature,\n    density,\n    Nt,\n    Δt,\n    seed,\n    gr_lastn\n)\n\ndefined at /home/runner/work/HappyMolecules.jl/HappyMolecules.jl/src/applications.jl:31.\n\n\n\n\n\n","category":"method"},{"location":"#HappyMolecules.jl","page":"Home","title":"HappyMolecules.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This repo is for tutorial purpose!","category":"page"},{"location":"","page":"Home","title":"Home","text":"Please check the Chapter 4 of the following book:","category":"page"},{"location":"","page":"Home","title":"Home","text":"(Image: )","category":"page"},{"location":"#Notebooks","page":"Home","title":"Notebooks","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A demo notebook","category":"page"}]
}