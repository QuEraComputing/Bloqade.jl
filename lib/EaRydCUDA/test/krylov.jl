using Test
using CUDA
using CUDA.CUSPARSE: CuSparseMatrixCSC,
    CuSparseMatrixCSR,
    AbstractCuSparseMatrix
using YaoArrayRegister
using YaoSubspaceArrayReg
using LinearAlgebra
using EaRydWaveforms
using EaRydKrylov
using EaRydExpr
using EaRydCUDA
using Adapt
CUDA.allowscalar(false)

@testset "KrylovEvolution" begin
    atoms = [(i*3.2, ) for i in 1:5]
    clocks=[0.0, 0.5, 0.8, 1.1, 1.5]
    wf = piecewise_constant(clocks=clocks, values=[0.0, 2.1, 2.1, 1.5, 0.0])
    h = rydberg_h(atoms; Ω=wf)
    reg = zero_state(length(atoms))

    prob = KrylovEvolution(reg, clocks, h)
    d_prob = adapt(CuArray, prob)
    emulate!(d_prob)
    emulate!(prob)
    @test Array(d_prob.reg.state) ≈ prob.reg.state
end
