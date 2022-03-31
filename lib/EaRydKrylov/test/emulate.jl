using Test
using EaRydExpr
using YaoArrayRegister
using EaRydWaveforms
using EaRydKrylov
using EaRydLattices
using EaRydExpr: Hamiltonian

atoms = generate_sites(ChainLattice(), 4, scale=3.2);
wf = piecewise_constant(clocks=[0.0, 0.5, 0.8, 1.1, 1.5], values=[0.0, 2.1, 2.1, 1.5, 0.0])
h = rydberg_h(atoms; Ω=wf)
reg = zero_state(length(atoms))
prob = KrylovEvolution(reg, diff(wf.f.clocks), h)
emulate!(prob)
