using Test
using RydbergEmulator
using LightGraphs: SimpleGraph, add_edge!
using SparseArrays
using OrderedCollections
using CUDA

@testset "simple graph hamiltonian subspace" begin
    subspace = Subspace(test_graph)
    @test collect(keys(subspace)) == sort(test_subspace_v)
    @test collect(values(subspace)) == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

    # test hamiltonian creation
    Ω = Float64[2.5, 3.4, 0.2, 1.7, 4.3]
    Δ = Float64[1.2, 3.4, 2.6, 0.2, 1.8]
    ϕ = Float64[0.5, 0.4, -0.2, -1.2, 10.2]

    target = create_test_hamiltonian(Δ, Ω, ϕ)
    h = XTerm(Ω, ϕ) + ZTerm(Δ)
    @test SparseMatrixCSC(h, subspace) ≈ target

    H = SparseMatrixCSC(h, subspace)
    @test H ≈ update_term!(copy(H), h, subspace)

    @testset "cuda" begin
        if CUDA.functional()
            using CUDA
            using CUDA.CUSPARSE
            dΩ = cu(Ω)
            dϕ = cu(ϕ)
            dΔ = cu(Δ)
            h = XTerm(dΩ, dϕ) + ZTerm(dΔ)
            dH = CuSparseMatrixCSR(H)
            ds = cu(subspace)
            update_term!(dH, h, ds)
            @test isapprox(SparseMatrixCSC(dH), H; rtol=1e-7)
        end
    end
end
