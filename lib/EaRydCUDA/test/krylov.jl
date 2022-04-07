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
    @show typeof(t)
    f(0.1) * t
end

using LinearAlgebra
lmul!(1.0, d_prob.hamiltonian.ts[1])
@which 1.0 * d_prob.hamiltonian.ts[2]

b = Broadcast.broadcasted(*, 1.0, prob.hamiltonian.ts[1])
@code_lowered Base.Broadcast.instantiate(b)

function LinearAlgebra.opnorm(h::StepHamiltonian, p=2)
    H = sum(zip(h.h.fs, h.h.ts)) do (f, t)
        f(h.t) * t
    end
    return opnorm(H, p)
end

space = Subspace(5, [0, 1, 4, 8])
zero_state(space)|>typeof
emulate!(prob)
