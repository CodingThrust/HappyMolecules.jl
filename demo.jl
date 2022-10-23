### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ fc7c89c4-52c2-11ed-0593-df19e379a6d9
using Pkg; Pkg.activate()

# ╔═╡ 12d6e3cf-b79d-4616-b842-a4d7479fadae
using Revise

# ╔═╡ a37aa6cd-6b93-4fa5-8345-eccf6e95c8e9
using HappyMolecules

# ╔═╡ 6bc04ad5-264f-4732-bc4f-48ab46d3435d
using GLMakie

# ╔═╡ adf2706f-e69e-41cc-a205-3be418501488
using PlutoUI

# ╔═╡ 55a5cba9-9fd6-4740-abff-6830383e7e23
using StaticArrays

# ╔═╡ 8fe06450-ee59-42e1-b991-803b9ea8eecc
using Makie.Colors

# ╔═╡ 6bfbbe53-2b3b-401f-8f5b-02527e12c993
natoms = 100

# ╔═╡ cd91b971-4820-428e-b452-bdfe9792e85d
x0 = rand(SVector{3, Float64}, natoms)

# ╔═╡ c5fb9f78-2e30-4147-a3ec-dfc3b3224313
box = PeriodicBox(SVector(1.0, 1.0, 1.0))

# ╔═╡ 409048a8-9e7f-447c-86c6-995e683d182b
rc2 = 0.1

# ╔═╡ 746f55a2-9f58-4793-ad0c-62d39c6b4551
temperature = 1.0

# ╔═╡ 870fd253-d6e9-4921-94cf-97e6395a2443
md = molecule_dynamics(x0, box, temperature, rc2, 0.00001)

# ╔═╡ e94c8411-a202-4ef1-a2a6-557c7c1883c8
step!(md);

# ╔═╡ 0b7b46eb-8b33-48fa-b4e7-60003003cd3a
let
	fig, ax, lineplot = scatter(getindex.(md.x, 1), getindex.(md.x,2))
	filename = tempname() * ".mp4"
	record(fig, filename, 1:100; framerate = 30) do i
		step!(md)
		scatter(getindex.(md.x, 1), getindex.(md.x,2))
	end
	LocalResource(filename)
end

# ╔═╡ 6e47a227-4895-4002-82e1-97292ada1903
md.x

# ╔═╡ Cell order:
# ╠═fc7c89c4-52c2-11ed-0593-df19e379a6d9
# ╠═12d6e3cf-b79d-4616-b842-a4d7479fadae
# ╠═a37aa6cd-6b93-4fa5-8345-eccf6e95c8e9
# ╠═6bc04ad5-264f-4732-bc4f-48ab46d3435d
# ╠═adf2706f-e69e-41cc-a205-3be418501488
# ╠═6bfbbe53-2b3b-401f-8f5b-02527e12c993
# ╠═55a5cba9-9fd6-4740-abff-6830383e7e23
# ╠═cd91b971-4820-428e-b452-bdfe9792e85d
# ╠═c5fb9f78-2e30-4147-a3ec-dfc3b3224313
# ╠═409048a8-9e7f-447c-86c6-995e683d182b
# ╠═746f55a2-9f58-4793-ad0c-62d39c6b4551
# ╠═870fd253-d6e9-4921-94cf-97e6395a2443
# ╠═e94c8411-a202-4ef1-a2a6-557c7c1883c8
# ╠═8fe06450-ee59-42e1-b991-803b9ea8eecc
# ╠═0b7b46eb-8b33-48fa-b4e7-60003003cd3a
# ╠═6e47a227-4895-4002-82e1-97292ada1903
