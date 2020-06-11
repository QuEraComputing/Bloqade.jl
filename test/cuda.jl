using RydbergEmulator
using CUDA
using LightGraphs
using BenchmarkTools

T = Float32
graph = unit_disk_graph(lattice_atoms(10, 0.8, "square"), 1.0)
ϕs = rand(Complex{T}, 5)
hs = SimpleRydberg.(ϕs)
ts = rand(length(hs))
# prepair a zero state
subspace_v = subspace(graph)
reg = RydbergEmulator.zero_state(10, subspace_v)
dreg = cu(reg)
qaoa = QAOA{nv(graph)}(cu(subspace_v), hs, ts)

dreg |> 

QAOA{nv(graph)}(subspace_v, hs, ts)
