using Test
using EaRydExpr

function to_matrix(h::Hamiltonian, t::Real)
    return sum(zip(h.fs, h.ts)) do (f, term)
        f(t) * term
    end
end

using YaoBlocks

EaRydExpr.lower_expr(ComplexF64, SumOfXPhase(5, 0.1, 0.2), fullspace)

atoms = [(0.1, 0.2), (0.2, 0.3), (0.3, 0.4), (0.1, 0.5), (0.2, 0.1)]
mat(RydInteract(atoms, 0.5))

h = RydInteract(atoms, 0.5) + SumOfX(5, 1.2) - SumOfN(5, sin)

lh = EaRydExpr.Hamiltonian(ComplexF64, h, fullspace)
lh.fs

lh.ts[1][2, 2]
using LinearAlgebra
out_state = rand(ComplexF64, 1<<5)
in_state = rand(ComplexF64, 1<<5)
mul!(out_state, lh, in_state)
