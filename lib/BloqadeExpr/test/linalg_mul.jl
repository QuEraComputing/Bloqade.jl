using Test
using BloqadeExpr
using YaoAPI
using LinearAlgebra
using Random
using SparseArrays
using SparseMatricesCSR
using LuxurySparse
using BloqadeExpr: Hamiltonian

atoms = [(1, 1), (1, 2), (1, 3), (1, 4), (1, 5)]
h = rydberg_h(atoms; Ω = 1.0, Δ = sin)

hlist = Hamiltonian(Float64, h)
@test length(hlist.fs) == 2
@test hlist.fs[end] === one

H = sum(zip(hlist.fs, hlist.ts)) do (f, t)
    return f(0.1) * t
end

C = zeros(ComplexF64, 1 << 5)
B = rand(ComplexF64, 1 << 5)
ot = zeros(ComplexF64, 1 << 5)
mul!(ot, hlist(0.1), B)
@test  ot ≈ H * B

# ThreadedMatrix unit tests
  # ThreadedMatrix with SparseMatrixCSC
  # Threaded Matrix with CSR and Diagonal matrix inside

@testset "ThreadedMatrix for Parallel Emulation mul" begin

    # these values are not mutated and can be re-used throughout the unit tests
    B = rand(ComplexF64, 10)
    α = 0.5 + 1.0im
    β = 2.1 + 3.3im

    @testset "Transposed SparseMatrixCSC Matrix (mul)" begin

        C_original = rand(ComplexF64, 10)
        C_copy = deepcopy(C_original)

        Asrc = sprand(ComplexF64, 10, 10, 0.1)
        A = transpose(Asrc)
        A_threaded = BloqadeExpr.ThreadedMatrix(transpose(Asrc))


        SparseArrays.mul!(C_original, A, B, α, β) # C -> A B α + C β 
        BloqadeExpr.mul!(C_copy, A_threaded, B, α, β)

        @test C_original == C_copy
    end

    @testset "CSR-format Matrix (mul)" begin

        C_original = rand(ComplexF64, 10)
        C_copy = deepcopy(C_original)

        A = sprand(ComplexF64, 10, 10, 0.1)
        A_threaded = A |> transpose |> SparseMatrixCSR |>  BloqadeExpr.ThreadedMatrix

        SparseArrays.mul!(C_original, transpose(A), B, α, β) # C -> A B α + C β 
        BloqadeExpr.mul!(C_copy, A_threaded, B, α, β)

        @test C_original == C_copy
    end

    @testset "Diagonal Matrix (mul)" begin

        C_original = rand(ComplexF64, 10)
        C_copy = deepcopy(C_original)

        A = Diagonal(rand(ComplexF64, 10))
        A_threaded = A |> BloqadeExpr.ThreadedMatrix

        SparseArrays.mul!(C_original, A, B, α, β) # C -> A B α + C β 
        BloqadeExpr.mul!(C_copy, A_threaded, B, α, β)

        @test C_original == C_copy
    end

    @testset "Permutation Matrix (mul)" begin

        C_original = rand(ComplexF64, 10)
        C_copy = deepcopy(C_original)

        A = PermMatrix(shuffle(1:10), rand(ComplexF64, 10))
        A_threaded = A |> BloqadeExpr.ThreadedMatrix
       
        mul!(C_original, A, B, α, β)
        BloqadeExpr.mul!(C_copy, A_threaded, B, α, β)

        @test C_original == C_copy
    end

end