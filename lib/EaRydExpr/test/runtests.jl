using Test
using EaRydExpr

function to_matrix(h::Hamiltonian, t::Real)
    return sum(zip(h.fs, h.ts)) do (f, term)
        f(t) * term
    end
end

using YaoBlocks

