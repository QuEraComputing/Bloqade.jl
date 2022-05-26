using Test
using Random
using YaoSubspaceArrayReg
using Yao
using SparseArrays
using LinearAlgebra

@testset "apply" begin
    space = Subspace(10, sort(randperm(1 << 10)[1:76] .- 1))
    r = SubspaceArrayReg(randn(ComplexF64, 76), space)
    for g in [
        put(10, 2 => X),
        put(10, 3 => Rx(0.4)),
        put(10, 2:6 => matblock(rand_unitary(32))),
        control(10, (3, -5), (2, 7) => matblock(rand_unitary(4))),
    ]
        M = SparseMatrixCSC(mat(g))
        r2 = apply(r, g)
        @test r2.state ≈ M[space, space] * r.state
    end
end

@testset "expect" begin
    space = Subspace(10, sort(randperm(1 << 10)[1:76] .- 1))
    r = SubspaceArrayReg(randn(ComplexF64, 76), space)
    @test r' * r ≈ norm(r.state)^2
    for g in [put(10, 2 => X), control(10, (3, -5), (2, 7) => matblock(rand_unitary(4)))]
        M = SparseMatrixCSC(mat(g))
        @test expect(g, r) ≈ r.state' * M[space, space] * r.state
    end
end
