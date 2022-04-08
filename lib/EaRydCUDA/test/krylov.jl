using CUDA
using CUDA.CUSPARSE: CuSparseMatrixCSC,
    CuSparseMatrixCSR,
    AbstractCuSparseMatrix
using YaoArrayRegister
using YaoSubspaceArrayReg
using EaRydWaveforms
using EaRydKrylov
using EaRydExpr
using EaRydCUDA

CUDA.allowscalar(false)

atoms = [(i*3.2, ) for i in 1:5]
clocks=[0.0, 0.5, 0.8, 1.1, 1.5]
wf = piecewise_constant(clocks=clocks, values=[0.0, 2.1, 2.1, 1.5, 0.0])
h = rydberg_h(atoms; Ω=wf)
reg = zero_state(length(atoms))
prob = KrylovEvolution(reg, clocks, h)
d_prob = cu(prob)
emulate!(d_prob)

H = sum(zip(d_prob.hamiltonian.fs, d_prob.hamiltonian.ts)) do (f, t)
    f(0.1) * t
end


st = statevec(cu(reg))
dA = d_prob.hamiltonian.ts[1]
y = similar(st)

using StructArrays

StructArray{ComplexF32}((real(st), imag(st)))