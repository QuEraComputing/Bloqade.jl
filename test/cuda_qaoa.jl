using RydbergEmulator
using CUDA
using Yao
using CUDA.CUSPARSE
using LightGraphs
using BenchmarkTools
using Test
using ExponentialUtilities
using LinearAlgebra
using RydbergEmulator: subspace

N = 25
graph = unit_disk_graph(lattice_atoms(N, 0.8, "square"), 1.0)
ϕs = rand(5)
hs = SimpleRydberg.(ϕs)
ts = rand(length(hs))
# prepair a zero state
subspace_v = subspace(graph)
reg = RydbergEmulator.zero_state(N, subspace_v)
dreg = cu(reg)
dqaoa = CuQAOA{nv(graph)}(subspace_v, hs, ts)
qaoa = QAOA{nv(graph)}(subspace_v, hs, ts)

t1 = @benchmark CUDA.@sync(begin
    st = vec(dreg.state)
    update_hamiltonian!(dqaoa.h_cache, N, dqaoa.device_subspace_v, 1.0, ϕs[1]);
    arnoldi!(dqaoa.Ks, dqaoa.h_cache, st; ishermitian=true)
    # @show dqaoa.Ks.m
    expv!(st, 1.0, dqaoa.Ks; dexpHe=dqaoa.dexpHe)
end) setup=(dreg=cu(RydbergEmulator.zero_state($N, $subspace_v)))

t2 = @benchmark begin
    st = vec(reg.state)
    update_hamiltonian!(qaoa.hamiltonian_cache, N, qaoa.subspace_v, 1.0, ϕs[1])
    arnoldi!(qaoa.Ks, qaoa.hamiltonian_cache, st; ishermitian=true)
    expv!(vec(reg.state), 1.0, qaoa.Ks)
end setup=(reg=RydbergEmulator.zero_state($N, $subspace_v))

minimum(t2).time / minimum(t1).time

# CUDA.@sync hm = update_hamiltonian!(dqaoa.h_cache, N, dqaoa.device_subspace_v, 1.0, ϕs[1]);

# update_hamiltonian!(qaoa.hamiltonian_cache, N, qaoa.subspace_v, 1.0, ϕs[1])

# lmul!(-im, dqaoa.h_cache.nzVal)
# r1 = Yao.apply!(copy(dreg), dqaoa)


# st = vec(dreg.state);
# dqaoa.Ks.beta = zero(dqaoa.Ks.beta)
# fill!(dqaoa.Ks.H, 0)
# update_hamiltonian!(dqaoa.h_cache, N, dqaoa.device_subspace_v, 1.0, ϕs[1]);
# arnoldi!(dqaoa.Ks, dqaoa.h_cache, st; ishermitian=true);
# expv!(st, 1.0, dqaoa.Ks; dexpHe=dqaoa.dexpHe);

# dqaoa.Ks

# @which expv!(st, 1.0, dqaoa.Ks; dexpHe=dqaoa.dexpHe)
# dqaoa = CuQAOA{nv(graph)}(subspace_v, hs, ts)

# dt1 = @benchmark CUDA.@sync update_hamiltonian!(dqaoa.h_cache, N, dqaoa.device_subspace_v, 1.0, ϕs[1])
# t2 = @benchmark update_hamiltonian!(qaoa.hamiltonian_cache, N, qaoa.subspace_v, 1.0, ϕs[1])

# dt2 = @benchmark CUDA.@sync arnoldi!(dqaoa.Ks, dqaoa.h_cache, st; ishermitian=true)
# t2 = @benchmark arnoldi!(qaoa.Ks, qaoa.hamiltonian_cache, vec(reg.state); ishermitian=true)

# dt3 = @benchmark CUDA.@sync expv!(st, 1.0, dqaoa.Ks; dexpHe=dqaoa.dexpHe)
# t3 = @benchmark expv!(vec(reg.state), 1.0, qaoa.Ks)

# (minimum(t1).time + minimum(t2).time + minimum(t3).time)/(minimum(dt1).time + minimum(dt2).time + minimum(dt3).time)

