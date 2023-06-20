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


@test tr(hlist(0.1)) ≈ tr(H)

# ThreadedMatrix unit tests
  # ThreadedMatrix with SparseMatrixCSC
  # Threaded Matrix with CSR and Diagonal matrix inside

@testset "ThreadedMatrix for Parallel Emulation (tr)" begin

    # these values are not mutated and can be re-used throughout the unit tests

    @testset "Transposed SparseMatrixCSC Matrix (tr)" begin

        #C_original = rand(ComplexF64, 10)
        #C_copy = deepcopy(C_original)

        A = sprand(ComplexF64, 20, 20, 0.3)
        A_dense = Matrix(A)
        A_threaded = BloqadeExpr.ThreadedMatrix(transpose(A))
        
        @test tr(A_dense) ≈ BloqadeExpr.tr(A_threaded)
        
    end

    @testset "CSR-format Matrix (tr)" begin

        
        A = sprand(ComplexF64, 20, 20, 0.3)
        A_dense = Matrix(A)
        A_threaded = A |> transpose |> SparseMatrixCSR |>  BloqadeExpr.ThreadedMatrix

        @test tr(transpose(A_dense)) ≈ BloqadeExpr.tr(A_threaded)

    end

    @testset "Diagonal Matrix (tr)" begin


        A = Diagonal(rand(ComplexF64, 20))
        A_dense = Matrix(A)
        A_threaded = A |> BloqadeExpr.ThreadedMatrix

        @test tr(A_dense) ≈ BloqadeExpr.tr(A_threaded)

    end

    @testset "Permutation Matrix (tr)" begin


        A = PermMatrix(shuffle(1:20), rand(ComplexF64, 20))
        A_dense = Matrix(A)
        A_threaded = A |> BloqadeExpr.ThreadedMatrix
       
        @test tr(A_dense) ≈ BloqadeExpr.tr(A_threaded)
        
    end

end