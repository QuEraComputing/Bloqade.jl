### A Pluto.jl notebook ###
# v0.11.11

using Markdown
using InteractiveUtils

# ╔═╡ b96793ec-eda7-11ea-237b-2daaf5170aeb
using Pkg; Pkg.status()

# ╔═╡ 895e5adc-e161-11ea-0dbf-cd1b5aeb5956
using RydbergEmulator, Plots, MISExperimentUtils, Yao

# ╔═╡ a23a135e-e3d9-11ea-3452-a3cbd8793251
md"""
# Rydberg Hamiltonian

```math
H_{ryd} = \sum_{i,j} \frac{C}{|r_i - r_j|^6} n_i n_j + \sum_i \Omega_i \cdot  (e^{iϕ_i}|0⟩⟨1| + e^{-iϕ_i}|1⟩⟨0|) + \sum_i Δ_i σ_z
```

## Time Independent Emulation

emulate the following matrix expression

```math
\exp(-iH_{ryd}t)|\psi⟩
```

this can be simulated efficiently (not in poly, still exponential) via projecting the ``H_{ryd}`` in Krylov subspace ``(H_{ryd}, |ψ⟩)``.

## Blockade Approximation

The interaction term will constraint the solution in a constraint space, which can be found by looking for maximal cliques of the (unit disk) interaction graph. This can reduce the simulation space.
"""

# ╔═╡ 6b52529e-e3ef-11ea-32eb-8322ff3746cb
md"# Use Case?"

# ╔═╡ a59a81c2-e167-11ea-3dfd-8d851e3116a3
md"# Define the graph"

# ╔═╡ 0d473c12-e162-11ea-1bd9-dfe48ce613fd
atoms = square_lattice(20, 0.8)

# ╔═╡ ed6880c6-e18f-11ea-2dba-f5319a439d01
@doc square_lattice

# ╔═╡ 0cffc54a-e166-11ea-16cd-196aac923257
graph = unit_disk_graph(atoms, 1.5)

# ╔═╡ 14e0bc60-e166-11ea-09b8-9d809ebcde21
viz_atoms(atoms, graph)

# ╔═╡ 1f88793c-e166-11ea-1a74-77d0c432a7c3
md"exact MIS = $(exact_solve_mis(graph))"

# ╔═╡ 33634d74-e166-11ea-0733-1de6e1898dc1
bit_string = let ts = rand(5)*10, ϕs = rand(5) * 2π
	final_state = qaoa_on_graph(graph, ts, ϕs)
end |> MISExperimentUtils.measure!

# ╔═╡ 5ff566a6-e166-11ea-398d-af937ed47243
viz_config(atoms, graph, bit_string)

# ╔═╡ 9a3aac26-e167-11ea-0c3b-1fbd2e0c9231
md"# Define register and subspace"

# ╔═╡ 935155ea-e167-11ea-3b11-fb542dafbb2d
@doc Subspace

# ╔═╡ a16b7c6c-e166-11ea-338c-7df3e1f3806c
s = Subspace(graph)

# ╔═╡ 89858392-e167-11ea-08ae-054857d5c56a
@doc RydbergEmulator.zero_state

# ╔═╡ b0878b66-e167-11ea-2f2a-f3fc1cc63e4c
r = RydbergEmulator.zero_state(length(atoms), s)

# ╔═╡ 2dbafc84-e169-11ea-263c-9966c6898d64
-mean_nv(r)

# ╔═╡ b48763d4-e168-11ea-3666-bf811d7d44e4
md"# Define Hamiltonian"

# ╔═╡ 71658b80-e190-11ea-069d-2b9c3b2d5e5b
@doc simple_rydberg

# ╔═╡ ac1b0854-e168-11ea-2390-fdd849ad4903
hs = [simple_rydberg(length(atoms), rand()) for _ in 1:5]

# ╔═╡ f8c8fa3a-e168-11ea-279c-39ece6567f76
md"# Start Simulation"

# ╔═╡ d58f6e20-e168-11ea-3e81-fd3415adf4bd
@doc emulate!

# ╔═╡ 03a1dbc0-e169-11ea-1238-87eaee8f8f80
emulate!(r, rand(5), hs)

# ╔═╡ 4d9778b4-e169-11ea-0bc3-978fcae7dc3b
md"# Check Loss"

# ╔═╡ 520cdf26-e169-11ea-11e2-7735563a81f7
@doc mean_nv

# ╔═╡ 0d24d3a0-e169-11ea-137d-e584fb46c018
-mean_nv(r)

# ╔═╡ 5821a72a-e169-11ea-1347-15aea7b628dd
@doc soft_misloss

# ╔═╡ 5e325952-e169-11ea-300a-1f7901f1087d
soft_misloss(r, 0.3)

# ╔═╡ 68d7c160-e169-11ea-2d06-373aa675d00e
md"# Sample from register"

# ╔═╡ f2864d3a-e190-11ea-1a21-e9f69e275688
@doc measure

# ╔═╡ 701762b6-e169-11ea-35e4-dbcd02c974d8
samples = measure(r; nshots=1000)

# ╔═╡ 86dcad60-e169-11ea-1ffb-d702c3b91c8d
-mean_nv(samples)

# ╔═╡ 8ce117a2-e169-11ea-362f-5f3df236e5ca
soft_misloss(samples, 0.3)

# ╔═╡ afa1ed52-e169-11ea-3ee6-a304679b5f51
md"# More Advanced usage"

# ╔═╡ 2172bbfa-e16a-11ea-0983-99aaf5a66111
md"## Rydberg Hamiltonian with Interact Term"

# ╔═╡ 3669e29a-e16a-11ea-1fb3-3fd96aa1b646
RydInteract(1.0, atoms)

# ╔═╡ 5c4a11ce-e16a-11ea-3c8a-031ae4410d5e
XTerm(length(atoms), 2.0)

# ╔═╡ 710ff7a4-e16a-11ea-1979-25b280e4a457
XTerm(length(atoms), 2.0, 3.0)

# ╔═╡ 78d40d84-e16a-11ea-19fd-f3e0829d0c3e
ZTerm(length(atoms), 2.0)

# ╔═╡ 374f917e-e191-11ea-3ca1-e16dd339caaa
XTerm(rand(5), rand(5))

# ╔═╡ 812b9c76-e16a-11ea-19f4-d9fae526126d
interact_hs = [RydInteract(1.0, atoms) + XTerm(length(atoms), 1.0, rand()) for _ in 1:5]

# ╔═╡ d390cf70-e16a-11ea-1925-53d8432d731e
rr = RydbergEmulator.zero_state(length(atoms), s)

# ╔═╡ f0317df0-e16a-11ea-0599-7d2ff48a6d60
emulate!(rr, rand(5), interact_hs)

# ╔═╡ c2a91496-e16b-11ea-00c5-8d8c5d231bc5
zx_hs = [XTerm(length(atoms), 1.0, rand()) + ZTerm(length(atoms), 2.0) for _ in 1:5]

# ╔═╡ d68c7bba-e16b-11ea-3d2a-35f23b1e40e2
emulate!(rr, rand(5), zx_hs)

# ╔═╡ 0b574ffe-e16d-11ea-203c-79caeb9bf94e
md"## Full space simulation"

# ╔═╡ 14d7a9d4-e16d-11ea-3ff1-e728745af1a3
fullspace_r = Yao.zero_state(length(atoms))

# ╔═╡ 21861724-e16d-11ea-22ec-c7d8cd4651ee
emulate!(fullspace_r, rand(5), interact_hs)

# ╔═╡ 4c96f56e-e16d-11ea-06c6-e5a23fc2a4a8
state(fullspace_r)

# ╔═╡ 81ebc258-e3ef-11ea-33f8-49490dc3512c
md"""
# Post You User Cases

you can open an issue to post a description of your use case
"""

# ╔═╡ e6f8c426-e3ef-11ea-0061-5950049e62ec
md"run this in XSTM"

# ╔═╡ Cell order:
# ╠═b96793ec-eda7-11ea-237b-2daaf5170aeb
# ╠═895e5adc-e161-11ea-0dbf-cd1b5aeb5956
# ╟─a23a135e-e3d9-11ea-3452-a3cbd8793251
# ╟─6b52529e-e3ef-11ea-32eb-8322ff3746cb
# ╟─a59a81c2-e167-11ea-3dfd-8d851e3116a3
# ╠═0d473c12-e162-11ea-1bd9-dfe48ce613fd
# ╠═ed6880c6-e18f-11ea-2dba-f5319a439d01
# ╠═0cffc54a-e166-11ea-16cd-196aac923257
# ╠═14e0bc60-e166-11ea-09b8-9d809ebcde21
# ╟─1f88793c-e166-11ea-1a74-77d0c432a7c3
# ╟─33634d74-e166-11ea-0733-1de6e1898dc1
# ╠═5ff566a6-e166-11ea-398d-af937ed47243
# ╟─9a3aac26-e167-11ea-0c3b-1fbd2e0c9231
# ╠═935155ea-e167-11ea-3b11-fb542dafbb2d
# ╠═a16b7c6c-e166-11ea-338c-7df3e1f3806c
# ╠═89858392-e167-11ea-08ae-054857d5c56a
# ╠═b0878b66-e167-11ea-2f2a-f3fc1cc63e4c
# ╠═2dbafc84-e169-11ea-263c-9966c6898d64
# ╟─b48763d4-e168-11ea-3666-bf811d7d44e4
# ╠═71658b80-e190-11ea-069d-2b9c3b2d5e5b
# ╠═ac1b0854-e168-11ea-2390-fdd849ad4903
# ╟─f8c8fa3a-e168-11ea-279c-39ece6567f76
# ╠═d58f6e20-e168-11ea-3e81-fd3415adf4bd
# ╠═03a1dbc0-e169-11ea-1238-87eaee8f8f80
# ╠═4d9778b4-e169-11ea-0bc3-978fcae7dc3b
# ╠═520cdf26-e169-11ea-11e2-7735563a81f7
# ╠═0d24d3a0-e169-11ea-137d-e584fb46c018
# ╠═5821a72a-e169-11ea-1347-15aea7b628dd
# ╠═5e325952-e169-11ea-300a-1f7901f1087d
# ╟─68d7c160-e169-11ea-2d06-373aa675d00e
# ╠═f2864d3a-e190-11ea-1a21-e9f69e275688
# ╠═701762b6-e169-11ea-35e4-dbcd02c974d8
# ╠═86dcad60-e169-11ea-1ffb-d702c3b91c8d
# ╠═8ce117a2-e169-11ea-362f-5f3df236e5ca
# ╟─afa1ed52-e169-11ea-3ee6-a304679b5f51
# ╟─2172bbfa-e16a-11ea-0983-99aaf5a66111
# ╠═3669e29a-e16a-11ea-1fb3-3fd96aa1b646
# ╠═5c4a11ce-e16a-11ea-3c8a-031ae4410d5e
# ╠═710ff7a4-e16a-11ea-1979-25b280e4a457
# ╠═78d40d84-e16a-11ea-19fd-f3e0829d0c3e
# ╠═374f917e-e191-11ea-3ca1-e16dd339caaa
# ╠═812b9c76-e16a-11ea-19f4-d9fae526126d
# ╠═d390cf70-e16a-11ea-1925-53d8432d731e
# ╠═f0317df0-e16a-11ea-0599-7d2ff48a6d60
# ╠═c2a91496-e16b-11ea-00c5-8d8c5d231bc5
# ╠═d68c7bba-e16b-11ea-3d2a-35f23b1e40e2
# ╟─0b574ffe-e16d-11ea-203c-79caeb9bf94e
# ╠═14d7a9d4-e16d-11ea-3ff1-e728745af1a3
# ╠═21861724-e16d-11ea-22ec-c7d8cd4651ee
# ╠═4c96f56e-e16d-11ea-06c6-e5a23fc2a4a8
# ╟─81ebc258-e3ef-11ea-33f8-49490dc3512c
# ╟─e6f8c426-e3ef-11ea-0061-5950049e62ec
