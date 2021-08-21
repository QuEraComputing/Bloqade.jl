using Test
using CUDA
using CuRydbergEmulator
using CUDA.CUSPARSE: CuSparseMatrixCSC, CuSparseMatrixCSR, AbstractCuSparseMatrix
using RydbergEmulator: AbstractTerm
using SparseArrays

CUDA.allowscalar(false)

@testset "update_term" begin
    atoms = square_lattice(5, 0.8)
    space = blockade_subspace(atoms, 1.5)
    h = RydInteract(atoms) + XTerm(length(atoms), 1.0) - NTerm(length(atoms), 1.2)

    H = SparseMatrixCSC{ComplexF32}(h, space)
    cuH = CuSparseMatrixCSR{ComplexF32}(H)
    update_term!(cuH, h, cu(space.subspace_v))
    @test SparseMatrixCSC(cuH) ≈ H

    H = SparseMatrixCSC{ComplexF32}(h)
    cuH = CuSparseMatrixCSR{ComplexF32}(H)
    update_term!(cuH, h)
    @test SparseMatrixCSC(cuH) ≈ H
end

# st = rand(ComplexF32, size(H, 1))
# @benchmark $H * $st
# using BenchmarkTools
# st = CUDA.rand(ComplexF32, size(cuH, 1))

# cuH = CuSparseMatrixCSR{ComplexF32}(H)
# @benchmark CUDA.@sync $cuH * $st

# cuH = CuSparseMatrixCSC{ComplexF32}(H)
# @benchmark CUDA.@sync $cuH * $st
