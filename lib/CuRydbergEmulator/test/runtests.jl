using Test
using CUDA
using Adapt
using ContinuousEmulator
using CuRydbergEmulator
using SparseArrays

using CUDA.CUSPARSE: CuSparseMatrixCSC, CuSparseMatrixCSR, AbstractCuSparseMatrix
using RydbergEmulator: AbstractTerm
using ContinuousEmulator: ShordingerEquation, update_dstate!

CUDA.allowscalar(false)

atoms = square_lattice(5, 0.8)
space = blockade_subspace(atoms, 1.5)
h = RydInteract(atoms) + XTerm(length(atoms), 1.0) - NTerm(length(atoms), 1.2)

@testset "update_term" begin
    H = SparseMatrixCSC{ComplexF32}(h, space)
    cuH = CuSparseMatrixCSR{ComplexF32}(H)
    update_term!(cuH, h, cu(space))
    @test SparseMatrixCSC(cuH) ≈ H

    H = SparseMatrixCSC{ComplexF32}(h)
    cuH = CuSparseMatrixCSR{ComplexF32}(H)
    update_term!(cuH, h)
    @test SparseMatrixCSC(cuH) ≈ H
end

@testset "update_dstate L=$L" for L in [10, 300]
    state = CUDA.rand(L, 2)
    dstate = CUDA.zeros(L, 2)

    update_dstate!(dstate, state, RealLayout())
    host_state = Array(state)
    host_dstate = Array(dstate)
    ref_state = -im * (host_state[:, 1] + im * host_state[:, 2])
    @test host_dstate[:, 1] ≈ real(ref_state)
    @test host_dstate[:, 2] ≈ imag(ref_state)
end


@testset "emulate h=$name" for (name, h) in [
    "x+z" => XTerm(5, 1.0) + ZTerm(5, sin),
    "rydberg" => rydberg_h(atoms, sin, nothing, cos),
]

    @testset "T=$T" for T in [ComplexF32, ComplexF64]
        ref_state = zero_state(T, length(atoms), space)
        dcomplex_r = cu(ref_state)
        dreal_r = cu(zero_state(T, length(atoms), space, RealLayout()))

        emulate!(ref_state, 1e-3, h)
        emulate!(dcomplex_r, 1e-3, h)
        emulate!(dreal_r, 1e-3, h)
        @test cpu(dreal_r) ≈ ref_state
        @test cpu(dcomplex_r) ≈ ref_state
    end
end
