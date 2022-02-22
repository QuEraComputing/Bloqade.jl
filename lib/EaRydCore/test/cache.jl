using Test
using EaRydCore
using SparseArrays
using LinearAlgebra
using EaRydCore: split_const_term

atoms = square_lattice(4, 0.8)

@testset "split_const_term $(nameof(typeof(space)))" for space in [FullSpace(), blockade_subspace(atoms)]
    for h in [
        rydberg_h(atoms; Δ=0.1, Ω=0.1),
        rydberg_h(atoms; Δ=0.1, Ω=sin),
        rydberg_h(atoms; Δ=cos, Ω=sin),
        rydberg_h(atoms; Δ=cos, Ω=[sin, sin, sin, sin]),
        rydberg_h(atoms; Δ=[cos, cos, cos, cos], Ω=[sin, sin, sin, sin]),
    ]

        H = SparseMatrixCSC{ComplexF64}(h(0.1), space)
        tc = split_const_term(ComplexF64, h, space)
        M = sum(zip(tc.fs, tc.hs)) do (f, h)
            if h isa AbstractBlock
                f(0.1) * mat(h)
            else
                f(0.1) * h
            end
        end

        @test M ≈ H
    end
end
