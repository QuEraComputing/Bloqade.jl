using Test, Random, RydbergEmulator
using Graphs, Yao, SparseArrays
using LinearAlgebra

@testset "apply" begin
    ss = independent_set_subspace(smallgraph(:petersen))
    r = RydbergReg(randn(ComplexF64, 76), ss)
    for g in [put(10, 2=>X), put(10, 3=>Rx(0.4)), put(10, 2:6=>matblock(rand_unitary(32))), control(10, (3,-5), (2,7)=>matblock(rand_unitary(4)))]
        M = SparseMatrixCSC(mat(g))
        r2 = apply(r, g)
        @test r2.state ≈ M[ss, ss] * r.state
    end
end

@testset "expect" begin
    ss = independent_set_subspace(smallgraph(:petersen))
    r = RydbergReg(randn(ComplexF64, 76), ss)
    @test r' * r ≈ norm(r.state)^2
    for g in [put(10, 2=>X), control(10, (3,-5), (2,7)=>matblock(rand_unitary(4)))]
        M = SparseMatrixCSC(mat(g))
        @test expect(g, r) ≈ r.state' * M[ss, ss] * r.state
    end
end