using Test
using RydbergEmulator
using LightGraphs: SimpleGraph, add_edge!
using SparseArrays
using OrderedCollections
using LuxurySparse
using BitBasis
using Printf
using Yao
using Yao.ConstGate: P0, P1

function rydinteract(atoms, C)
    n = length(atoms)
    terms = []
    for i in 1:n, j in 1:n
        if i != j
            push!(terms, C/(2 * RydbergEmulator.distance(atoms[i], atoms[j])^6) * kron(n, i=>P1, j=>P1))
        end
    end
    return sum(terms)
end

function test_print(f, h)
    ans = sprint((io, x)->show(io, MIME("text/plain"), x), h; context=:color=>true)
    @test ans == f()
end

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

    # @testset "cuda" begin
    #     if CUDA.functional()
    #         using CUDA
    #         using CUDA.CUSPARSE
    #         dΩ = cu(Ω)
    #         dϕ = cu(ϕ)
    #         dΔ = cu(Δ)
    #         h = XTerm(dΩ, dϕ) + ZTerm(dΔ)
    #         H = SparseMatrixCSC(h, subspace)
    #         dH = CuSparseMatrixCSR(H)
    #         ds = cu(subspace)
    #         update_term!(dH, cu(h), ds)
    #         @test isapprox(SparseMatrixCSC(dH), H; rtol=1e-7)
    #     end
    # end
end

@testset "X term" begin
    H = mat(sum([2.0 * kron(5, k=>Yao.X) for k in 1:5]))
    h = XTerm(5, 2.0)
    @test SparseMatrixCSC(h) ≈ H
    @test update_term!(SparseMatrixCSC(h), h) ≈ H
    test_print(h) do
        """
        XTerm
         ∑(n=1:5) \e[32m2.00\e[39m\e[94m σ^x\e[39m"""
    end

    Ωs = rand(5)
    h = XTerm(Ωs)
    H = mat(sum([Ω * kron(5, k=>Yao.X) for (k, Ω) in enumerate(Ωs)]))
    @test SparseMatrixCSC(h) ≈ H
    @test update_term!(SparseMatrixCSC(h), h) ≈ H
    @test eltype(h) == Float64
    test_print(h) do
        """
        XTerm
         \e[32m$(@sprintf("%.2f", Ωs[1]))\e[39m\e[94m σ^x\e[39m +
         \e[32m$(@sprintf("%.2f", Ωs[2]))\e[39m\e[94m σ^x\e[39m +
         \e[32m$(@sprintf("%.2f", Ωs[3]))\e[39m\e[94m σ^x\e[39m +
         \e[32m$(@sprintf("%.2f", Ωs[4]))\e[39m\e[94m σ^x\e[39m +
         \e[32m$(@sprintf("%.2f", Ωs[5]))\e[39m\e[94m σ^x\e[39m"""
    end

    Ωs = rand(5)
    ϕs = rand(5)
    h = XTerm(Ωs, ϕs)
    test_print(h) do
        """
        XTerm
         \e[32m$(@sprintf("%.2f", Ωs[1]))\e[39m (e^{\e[32m$(@sprintf("%.2f", ϕs[1]))\e[39mi}\e[94m|0)⟨1|\e[39m + e^{-\e[32m$(@sprintf("%.2f", ϕs[1]))\e[39mi}\e[94m|1⟩⟨0|\e[39m) +
         \e[32m$(@sprintf("%.2f", Ωs[2]))\e[39m (e^{\e[32m$(@sprintf("%.2f", ϕs[2]))\e[39mi}\e[94m|0)⟨1|\e[39m + e^{-\e[32m$(@sprintf("%.2f", ϕs[2]))\e[39mi}\e[94m|1⟩⟨0|\e[39m) +
         \e[32m$(@sprintf("%.2f", Ωs[3]))\e[39m (e^{\e[32m$(@sprintf("%.2f", ϕs[3]))\e[39mi}\e[94m|0)⟨1|\e[39m + e^{-\e[32m$(@sprintf("%.2f", ϕs[3]))\e[39mi}\e[94m|1⟩⟨0|\e[39m) +
         \e[32m$(@sprintf("%.2f", Ωs[4]))\e[39m (e^{\e[32m$(@sprintf("%.2f", ϕs[4]))\e[39mi}\e[94m|0)⟨1|\e[39m + e^{-\e[32m$(@sprintf("%.2f", ϕs[4]))\e[39mi}\e[94m|1⟩⟨0|\e[39m) +
         \e[32m$(@sprintf("%.2f", Ωs[5]))\e[39m (e^{\e[32m$(@sprintf("%.2f", ϕs[5]))\e[39mi}\e[94m|0)⟨1|\e[39m + e^{-\e[32m$(@sprintf("%.2f", ϕs[5]))\e[39mi}\e[94m|1⟩⟨0|\e[39m)"""
    end
end

@testset "Z term" begin
    H = mat(sum([2.0 * kron(5, k=>Yao.Z) for k in 1:5]))
    h = ZTerm(5, 2.0)
    @test SparseMatrixCSC(h) ≈ H
    @test update_term!(SparseMatrixCSC(h), h) ≈ H

    test_print(h) do
        """
        ZTerm
         ∑(n=1:5) \e[32m2.00\e[39m\e[94m σ^z\e[39m"""
    end

    Δs = rand(5)
    h = ZTerm(Δs)
    H = mat(sum([Δs[k] * kron(5, k=>Yao.Z) for k in 1:5]))
    @test SparseMatrixCSC(h) ≈ H
    @test update_term!(SparseMatrixCSC(h), h) ≈ H
    @test eltype(h) == Float64
    test_print(h) do
        """
        ZTerm
         \e[32m$(@sprintf("%.2f", Δs[1]))\e[39m\e[94m σ^z\e[39m +
         \e[32m$(@sprintf("%.2f", Δs[2]))\e[39m\e[94m σ^z\e[39m +
         \e[32m$(@sprintf("%.2f", Δs[3]))\e[39m\e[94m σ^z\e[39m +
         \e[32m$(@sprintf("%.2f", Δs[4]))\e[39m\e[94m σ^z\e[39m +
         \e[32m$(@sprintf("%.2f", Δs[5]))\e[39m\e[94m σ^z\e[39m"""
    end
end

@testset "rydberg interact term" begin
    atoms = RydbergEmulator.square_lattice(4, 0.8)
    H = mat(rydinteract(atoms, 2.0))
    h = RydInteract(atoms, 2.0)
    @test SparseMatrixCSC(h) ≈ H
    @test update_term!(SparseMatrixCSC(h), h) ≈ H
    @test eltype(h) == Float64
end

@testset "composite term" begin
    atoms = RydbergEmulator.square_lattice(5, 0.8)
    h = XTerm(5, 2.0) + RydInteract(atoms, 2.0) + ZTerm(5, 1.0)
    H = rydinteract(atoms, 2.0) +
        sum([1.0 * kron(5, k=>Yao.Z) for k in 1:5]) +
        sum([2.0 * kron(5, k=>Yao.X) for k in 1:5])
    H = mat(H)
    @test SparseMatrixCSC(h) ≈ H
    @test update_term!(SparseMatrixCSC(h), h) ≈ H
end

@testset "XTerm subspace" begin
    Ω = Float64[2.5, 3.4, 0.2, 1.7, 4.3]
    ϕ = Float64[0.5, 0.4, -0.2, -1.2, 10.2]
    h = XTerm(Ω, ϕ)

    subspace = Subspace(test_graph)
    H = SparseMatrixCSC(h, subspace)
    @test update_term!(copy(H), h, subspace) ≈ H
end

@testset "redberg interact term subspace" begin
    atoms = RydbergEmulator.square_lattice(4, 0.8)
    graph = unit_disk_graph(atoms, 1.5)
    s = Subspace(graph)
    H = mat(rydinteract(atoms, 2.0));
    h = RydInteract(atoms, 2.0)

    @test SparseMatrixCSC(h)[s.subspace_v.+1, s.subspace_v.+1] ≈ SparseMatrixCSC(h, s)
    @test update_term!(SparseMatrixCSC(h, s), h, s) ≈ SparseMatrixCSC(h, s)
end

@testset "utils" begin
    @test RydbergEmulator.to_tuple((1, 2, 3)) == (1, 2, 3)
    @test RydbergEmulator.to_tuple([1, 2, 3]) == (1, 2, 3)
end
