using BloqadeODE: BloqadeSolver, DP8Solver, SolverIterator, integrate!
using BloqadeExpr: rydberg_h, nqudits
using Bloqade
using YaoArrayRegister


h = rydberg_h([(0,0), (0,6)], Ω=t -> sin(t))
reg = zero_state(nqudits(h))


integrate!(BloqadeSolver(reg, 0, h), collect(0:0.01:1.0)) do t, state
    [rydberg_density(state, i) for i in 1:nqudits(h)]
end
