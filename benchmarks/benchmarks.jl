using BenchmarkTools
using CUDA
using CUDA.CUSPARSE
using RydbergEmulator
using LightGraphs
using SparseArrays
using Random

const SUITE = BenchmarkGroup()
update_suite = SUITE["update_term!"] = BenchmarkGroup()
emulate_suite = SUITE["emulate!"] = BenchmarkGroup()

Random.seed!(42)
n = 20
atoms = RydbergEmulator.square_lattice(n, 0.8)
g = unit_disk_graph(atoms)
s = Subspace(g)
h = simple_rydberg(n, 2.0)
H = SparseMatrixCSC(h)

update_suite["cpu"] = @benchmarkable update_term!($H, $h)

if CUDA.functional()
    dH = CuSparseMatrixCSR(H)
    update_suite["cuda"] = @benchmarkable CUDA.@sync update_term!($dH, $h)
end

Random.seed!(42)
ts = rand(5)
hs = simple_rydberg.(n, rand(5))
cache = EmulatorCache(ts, hs, s);
emulate_suite["simple_rydberg cpu"] = @benchmarkable emulate!(r, $ts, $hs, $cache) setup=(r = RydbergEmulator.zero_state($n, $s))

if CUDA.functional()
    dhs = cu.(hs)
    dcache = cu(cache)
    emulate_suite["simple_rydberg cuda"] = @benchmarkable CUDA.@sync(emulate!(r, $ts, $dhs, $dcache)) setup=(r = cu(RydbergEmulator.zero_state($n, $s)))
end
