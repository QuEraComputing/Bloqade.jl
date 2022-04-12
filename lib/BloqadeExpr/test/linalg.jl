using Test
using BloqadeExpr
using YaoAPI
using LinearAlgebra
using BloqadeExpr: Hamiltonian

atoms = [(1, 1), (1, 2), (1, 3), (1, 4), (1, 5)]
h = rydberg_h(atoms; Ω=1.0, Δ=sin)

hlist = Hamiltonian(Float64, h)
@test length(hlist.fs) == 2
@test hlist.fs[end] === one

H = sum(zip(hlist.fs, hlist.ts)) do (f,t)
    f(0.1) * t
end

C = zeros(ComplexF64, 1<<5)
B = rand(ComplexF64, 1<<5)
@test mul!(zeros(ComplexF64, 1<<5), hlist(0.1), B) ≈ H * B
