using EaRydCore
using SparseArrays
using EaRydCore: update_term!, add_term!

atoms = square_lattice(20, 0.8)
h = rydberg_h(atoms; Δ=0.1, Ω=0.2)
using BenchmarkTools
H = SparseMatrixCSC(h)
@benchmark update_term!(H, h)
@benchmark update_term!(H, h[1])
@benchmark update_term!(H, h[2])
@benchmark update_term!(H, h[3])
