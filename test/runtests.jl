using Test
using CUDA
using Adapt
using ContinuousEmulator
using CuRydbergEmulator
using SparseArrays

using CUDA.CUSPARSE: CuSparseMatrixCSC, CuSparseMatrixCSR, AbstractCuSparseMatrix
using RydbergEmulator: AbstractTerm
using ContinuousEmulator: ShordingerEquation


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

@testset "cu(ShordingerEquation)" begin
    eq = ShordingerEquation(Float32, h, space)
    @test eltype(eq.cache.state) == ComplexF32
    @test eltype(eq.cache.hamiltonian) == Float32

    # hamiltonian should be complex since cuSPARSE
    # doesn't support complex-real mv routine
    deq = cu(eq)
    @test eltype(deq.cache.state) == ComplexF32
    @test eltype(deq.cache.hamiltonian) == ComplexF32

    # we don't convert to Float32 automatically like cu(x)
    # we still want explicit control of numerical precision
    eq = ShordingerEquation(h, space)
    deq = cu(eq)
    @test eltype(deq.cache.state) == ComplexF64
    @test eltype(deq.cache.hamiltonian) == ComplexF64
end

@testset "ode solver" begin
    r = zero_state(ComplexF32, length(atoms), space)
    dr = cu(r)

    emulate!(r, 1e-3, h)
    emulate!(dr, 1e-3, h)

    @test cpu(dr) ≈ r

    r = zero_state(ComplexF64, length(atoms), space)
    dr = cu(r)

    emulate!(r, 1e-3, h)
    emulate!(dr, 1e-3, h)

    @test cpu(dr) ≈ r
end

# # st = rand(ComplexF32, size(H, 1))
# # @benchmark $H * $st
# # using BenchmarkTools
# # st = CUDA.rand(ComplexF32, size(cuH, 1))

# # cuH = CuSparseMatrixCSR{ComplexF32}(H)
# # @benchmark CUDA.@sync $cuH * $st

# # cuH = CuSparseMatrixCSC{ComplexF32}(H)
# # @benchmark CUDA.@sync $cuH * $st

# using BenchmarkTools
# 7 * 7 * 0.8
# atoms = square_lattice(40, 0.8)
# graph = unit_disk_graph(atoms, 1.5)
# space = blockade_subspace(graph)
# r = zero_state(ComplexF32, 40, space)

# h = RydInteract(atoms) + XTerm(length(atoms), 1.0) - NTerm(length(atoms), 1.2)
# r = zero_state(ComplexF32, length(atoms), space)
# dr = cu(r)

# SparseMatrixCSC(h, space)

# emulate!(r, 1e-3, h)
# @time emulate!(dr, 1e-3, h; progress=true, progress_steps=1)

# @benchmark emulate!($r, 1e-3, h)
# @benchmark CUDA.@sync emulate!($dr, 1e-3, h)

# using ContinuousEmulator
# ContinuousEvolution

# using CuRydbergEmulator
# using RydbergEmulator: PrecisionAdaptor
# using Adapt
# using CUDA
# x = CUDA.rand(10)
# adapt(PrecisionAdaptor(Float64), x)

# atoms = square_lattice(10, 0.8)
# space = blockade_subspace(atoms, 1.5)
# r = zero_state(10, space)
# dr = cu(r)
# h = rydberg_h(atoms, sin, nothing, 2.0)
# prob = ContinuousEvolution(dr, 2.0, h)
# emulate!(prob)
