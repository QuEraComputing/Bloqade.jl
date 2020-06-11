using RydbergEmulator
using CUDA
using LightGraphs
using BenchmarkTools

T = Float32
graph = unit_disk_graph(lattice_atoms(10, 0.8, "square"), 1.0)
ϕs = rand(Complex{T}, 5)
hs = SimpleRydberg.(ϕs)
# prepair a zero state
subspace_v = cu(subspace(graph))
st = CUDA.zeros(ComplexF32, length(subspace_v)); st[1] = 1
cureg = RydbergReg{nv(graph)}(st, subspace_v)

cureg |> QAOA{nv(graph)}(subspace_v, hs, ts)

QAOA{nv(graph)}(subspace_v, hs, ts)