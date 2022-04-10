using CUDA
using CUDA.CUSPARSE: CuSparseMatrixCSC,
    CuSparseMatrixCSR,
    AbstractCuSparseMatrix
using YaoArrayRegister
using ComplexArrays
using YaoSubspaceArrayReg
using EaRydWaveforms
using EaRydKrylov
using EaRydExpr
using EaRydCUDA
using Adapt
using BenchmarkTools
CUDA.allowscalar(false)

n = 10
atoms = [(i*3.2, ) for i in 1:n]
clocks=[0.0, 0.5, 0.8, 1.1, 1.5]
wf = piecewise_constant(clocks=clocks, values=[0.0, 2.1, 2.1, 1.5, 0.0])
h = rydberg_h(atoms; Ω=wf)
reg = zero_state(length(atoms))
reg = adapt(ComplexArray, reg)

prob = KrylovEvolution(reg, clocks, h)
d_prob = adapt(CuArray, prob)

@time emulate!(d_prob)
@time emulate!(prob)

Array(d_prob.reg.state.storage) ≈ prob.reg.state.storage

@profview emulate!(d_prob)

@benchmark emulate!(d_prob)
@benchmark emulate!(prob)


x = adapt(ComplexArray, rand(ComplexF64, 1<<n))
A = prob.hamiltonian(0.1)
y = similar(x)
d_prob = adapt(CuArray, prob)
dx = adapt(CuArray, x)
dA = d_prob.hamiltonian(0.1)
dy = similar(dx)

H = EaRydExpr.to_matrix(A)


using EaRydKrylov: expmv!
@time expmv!(dy, 0.1im, dA, dx)
@time expmv!( y, 0.1im,  A,  x)
