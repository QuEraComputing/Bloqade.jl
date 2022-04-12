using Test
using BloqadeKrylov
using SparseArrays

@testset "expmv" begin
    A = rand(10, 10)
    st = rand(10)
    @test BloqadeKrylov.expmv(-0.1im, A, st) ≈ exp(-0.1im * A) * st

    A = sprand(ComplexF64, 10, 10, 0.8)
    st = rand(ComplexF64, 10)
    @test BloqadeKrylov.expmv(-0.1im, A, st) ≈ exp(-0.1im * Matrix(A)) * st
end
