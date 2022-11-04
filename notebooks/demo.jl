### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ fc7c89c4-52c2-11ed-0593-df19e379a6d9
using Pkg; Pkg.activate()

# ╔═╡ 12d6e3cf-b79d-4616-b842-a4d7479fadae
using Revise

# ╔═╡ a37aa6cd-6b93-4fa5-8345-eccf6e95c8e9
using HappyMolecules

# ╔═╡ 6bc04ad5-264f-4732-bc4f-48ab46d3435d
using CairoMakie

# ╔═╡ adf2706f-e69e-41cc-a205-3be418501488
using PlutoUI

# ╔═╡ 55a5cba9-9fd6-4740-abff-6830383e7e23
using StaticArrays

# ╔═╡ 6bfbbe53-2b3b-401f-8f5b-02527e12c993
natoms = 100

# ╔═╡ c5fb9f78-2e30-4147-a3ec-dfc3b3224313
box = PeriodicBox(20.0, 20.0)

# ╔═╡ 0801e1e9-2468-48a2-aaae-5ed22fdd4783
x0 = uniform_locations(box, natoms)

# ╔═╡ 409048a8-9e7f-447c-86c6-995e683d182b
rc2 = 10.0

# ╔═╡ 746f55a2-9f58-4793-ad0c-62d39c6b4551
temperature = 1.0

# ╔═╡ 870fd253-d6e9-4921-94cf-97e6395a2443
md = molecule_dynamics(x0, box, temperature, rc2, 0.005)

# ╔═╡ e94c8411-a202-4ef1-a2a6-557c7c1883c8
step!(md);

# ╔═╡ 83abee9d-3aa4-4d39-984d-2b4561bc56d8
md"""
 $(@bind generate_video CheckBox()) generate clip
"""

# ╔═╡ 414610d2-c57e-440a-862d-7be78b9b4d37
filename = tempname() * ".mp4"

# ╔═╡ 0b7b46eb-8b33-48fa-b4e7-60003003cd3a
if generate_video
	points = Observable([Point2f(x...,) for x in md.x])
	directions = Observable([Point2f(x/10...,) for x in md.field])
	fig, ax = scatter(points)
	arrows!(ax, points, directions)
	limits!(ax, 0, box.dimensions[1], 0, box.dimensions[2])
	record(fig, filename, 1:500; framerate = 30, sleep=true) do i
		for j=1:10
			step!(md)
		end
		points[] = [Point2f(mod.(x, box.dimensions)...,) for x in md.x]
		directions[] = [Point2f(x/10...,) for x in md.field]
	end
	LocalResource(filename)
end

# ╔═╡ 1d31d9dd-42f5-49e3-a22c-70267f5442ce


# ╔═╡ Cell order:
# ╠═fc7c89c4-52c2-11ed-0593-df19e379a6d9
# ╠═12d6e3cf-b79d-4616-b842-a4d7479fadae
# ╠═a37aa6cd-6b93-4fa5-8345-eccf6e95c8e9
# ╠═6bc04ad5-264f-4732-bc4f-48ab46d3435d
# ╠═adf2706f-e69e-41cc-a205-3be418501488
# ╠═6bfbbe53-2b3b-401f-8f5b-02527e12c993
# ╠═55a5cba9-9fd6-4740-abff-6830383e7e23
# ╠═c5fb9f78-2e30-4147-a3ec-dfc3b3224313
# ╠═0801e1e9-2468-48a2-aaae-5ed22fdd4783
# ╠═409048a8-9e7f-447c-86c6-995e683d182b
# ╠═746f55a2-9f58-4793-ad0c-62d39c6b4551
# ╠═870fd253-d6e9-4921-94cf-97e6395a2443
# ╠═e94c8411-a202-4ef1-a2a6-557c7c1883c8
# ╟─83abee9d-3aa4-4d39-984d-2b4561bc56d8
# ╠═414610d2-c57e-440a-862d-7be78b9b4d37
# ╠═0b7b46eb-8b33-48fa-b4e7-60003003cd3a
# ╠═1d31d9dd-42f5-49e3-a22c-70267f5442ce
