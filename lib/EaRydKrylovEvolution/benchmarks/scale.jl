using BenchmarkTools
using EaRydKrylovEvolution
using SparseArrays
using Graphs
using Random
using CUDA
using CUDA.CUSPARSE
using UnicodePlots

N = 20
atoms = EaRydKrylovEvolution.square_lattice(N, 0.8)
h = XTerm(N, 0.4) + ZTerm(N, 0.3)
H = SparseMatrixCSC{ComplexF32}(h)
dH = CuSparseMatrixCSR(H)
dh = cu(h)

@benchmark update_term!($H, $h)
# BenchmarkTools.Trial: 
#   memory estimate:  0 bytes
#   allocs estimate:  0
#   --------------
#   minimum time:     123.915 ms (0.00% GC)
#   median time:      126.799 ms (0.00% GC)
#   mean time:        127.089 ms (0.00% GC)
#   maximum time:     132.306 ms (0.00% GC)
#   --------------
#   samples:          40
#   evals/sample:     1

@benchmark update_term!($dH, $dh)
# BenchmarkTools.Trial: 
#   memory estimate:  2.02 KiB
#   allocs estimate:  64
#   --------------
#   minimum time:     3.051 ms (0.00% GC)
#   median time:      3.161 ms (0.00% GC)
#   mean time:        3.215 ms (0.00% GC)
#   maximum time:     8.143 ms (0.00% GC)
#   --------------
#   samples:          1555
#   evals/sample:     1

update_term!(dH, dh)

N = 10
h = RydInteract(atoms, 2.0) + XTerm(N, 0.4) + ZTerm(N, 0.3)
H = SparseMatrixCSC{ComplexF32}(h);
dH = CuSparseMatrixCSR(H);

dh = cu(h)

@benchmark update_term!($H, $h)
@benchmark update_term!($dH, $dh)
update_term!(dH, dh)

CUDA.allowscalar(false)
CuSparseMatrixCSR(H);

function runbenchmarks(nd)
    cpu_benchmarks = []
    gpu_benchmarks = []
    for n in nd
        h = RydInteract(atoms, 2.0) + XTerm(N, 0.4) + ZTerm(N, 0.3)
        H = SparseMatrixCSC{ComplexF32}(h)
        dH = CuSparseMatrixCSR(H)
        dh = cu(h)

        cpu = @benchmarkable(update_term!($H, $h))
        gpu = @benchmarkable(update_term!($dH, $dh))

        push!(cpu_benchmarks, cpu)
        push!(gpu_benchmarks, gpu)
    end

    cpu_results = map(run, cpu_benchmarks)
    gpu_results = map(run, gpu_benchmarks)

    p_cpu = lineplot(nd, map(x->(minimum(x).time), cpu_results), title="cpu scaling")
    p_gpu = lineplot(nd, map(x->(minimum(x).time), gpu_results), title="gpu scaling")

    repr(p_cpu)
    repr(p_gpu)
    return cpu_results, gpu_results
end

cpu_results, gpu_results = runbenchmarks(8:22)

map(x->minimum(x).time, cpu_results)


subspace = Subspace(test_graph)
Ω = Float64[2.5, 3.4, 0.2, 1.7, 4.3]
Δ = Float64[1.2, 3.4, 2.6, 0.2, 1.8]
ϕ = Float64[0.5, 0.4, -0.2, -1.2, 10.2]
dΩ = cu(Ω)
dϕ = cu(ϕ)
dΔ = cu(Δ)
h = XTerm(dΩ, dϕ)# + ZTerm(dΔ)
H = SparseMatrixCSC(h)
dH = CuSparseMatrixCSR(H)
ds = cu(subspace)
update_term!(dH, h, ds)
