using CUDA
using CUDA.CUSPARSE: CuSparseMatrixCSC, CuSparseMatrixCSR, AbstractCuSparseMatrix
using YaoArrayRegister
using ComplexArrays
using YaoSubspaceArrayReg
using BloqadeWaveforms
using BloqadeKrylov
using BloqadeExpr
using BloqadeCUDA
using Adapt
using BenchmarkTools
CUDA.allowscalar(false)

function benchmark_suite(n)
    atoms = [(i * 3.2,) for i in 1:n]
    clocks = [0.0, 0.5, 0.8, 1.1, 1.5]
    wf = piecewise_constant(clocks = clocks, values = [0.0, 2.1, 2.1, 1.5, 0.0])
    h = rydberg_h(atoms; Ω = wf)
    reg = zero_state(length(atoms))
    prob = KrylovEvolution(reg, clocks, h)
    d_prob = adapt(CuArray, prob)
    @info "benchmarking cpu..."
    cpu = @elapsed emulate!(prob)
    @info "benchmarking cuda..."
    cuda = @elapsed CUDA.@sync emulate!(d_prob)
    return cpu, cuda
end

report = (cpu = Float64[], cuda = Float64[])
for n in 5:16
    @info "benchmarking..." n
    cpu, cuda = benchmark_suite(n)
    push!(report.cpu, cpu)
    push!(report.cuda, cuda)
end

report.cpu ./ report.cuda
