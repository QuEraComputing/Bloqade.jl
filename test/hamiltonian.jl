using Test
using RydbergEmulator
using LightGraphs
using ExponentialUtilities
using LinearAlgebra
using RydbergEmulator: subspace

function naive_qaoa(st, g, hs, ts)
    for (h, t) in zip(hs, ts)
        st = expv(t, -im * to_matrix(g, h.Ω, h.ϕ), st)
    end
    return st
end


@testset "hamiltonian" begin
    g = SimpleGraph(5)
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)
    add_edge!(g, 2, 4)
    add_edge!(g, 2, 5)
    add_edge!(g, 3, 4)
    add_edge!(g, 4, 5)

    # test maximal_cliques
    cg = complement(g)
    mis = maximal_cliques(cg)
    @test sort(mis) == [[1,3,5],[1,4],[2]]

    # test subspace
    set_subspace = [0,1,2,4,8,16,5,9,17,20,21]
    @test sort(subspace(5,mis)) == sort(set_subspace)

    # test to_matrix (hamiltonian creation)
    Ω = Float64[2.5, 3.4, 0.2, 1.7, 4.3]
    Δ = Float64[1.2, 3.4, 2.6, 0.2, 1.8]
    ϕ = Float64[0.5, 0.4, -0.2, -1.2, 10.2]

    # Ω = randn(Float64,5,1)
    # Δ = randn(Float64,5,1)
    # ϕ = randn(Float64,5,1)

    matrix_hamiltonian = zeros(ComplexF64,11,11)
    # diagonal sz terms
    matrix_hamiltonian[1,1] = Δ[1] + Δ[2] + Δ[3] + Δ[4] + Δ[5]
    matrix_hamiltonian[2,2] = - Δ[1] + Δ[2] + Δ[3] + Δ[4] + Δ[5]
    matrix_hamiltonian[3,3] = Δ[1] - Δ[2] + Δ[3] + Δ[4] + Δ[5]
    matrix_hamiltonian[4,4] = Δ[1] + Δ[2] - Δ[3] + Δ[4] + Δ[5]
    matrix_hamiltonian[5,5] = - Δ[1] + Δ[2] - Δ[3] + Δ[4] + Δ[5]
    matrix_hamiltonian[6,6] = Δ[1] + Δ[2] + Δ[3] - Δ[4] + Δ[5]
    matrix_hamiltonian[7,7] = - Δ[1] + Δ[2] + Δ[3] - Δ[4] + Δ[5]
    matrix_hamiltonian[8,8] = + Δ[1] + Δ[2] + Δ[3] + Δ[4] - Δ[5]
    matrix_hamiltonian[9,9] = - Δ[1] + Δ[2] + Δ[3] + Δ[4] - Δ[5]
    matrix_hamiltonian[10,10] = + Δ[1] + Δ[2] - Δ[3] + Δ[4] - Δ[5]
    matrix_hamiltonian[11,11] = - Δ[1] + Δ[2] - Δ[3] + Δ[4] - Δ[5]

    # spin flip terms
    matrix_hamiltonian[1,2] = Ω[1] * exp(im * ϕ[1])
    matrix_hamiltonian[2,1] = Ω[1] * exp(-im * ϕ[1])

    matrix_hamiltonian[1,3] = Ω[2] * exp(im * ϕ[2])
    matrix_hamiltonian[3,1] = Ω[2] * exp(-im * ϕ[2])

    matrix_hamiltonian[1,4] = Ω[3] * exp(im * ϕ[3])
    matrix_hamiltonian[4,1] = Ω[3] * exp(-im * ϕ[3])

    matrix_hamiltonian[1,6] = Ω[4] * exp(im * ϕ[4])
    matrix_hamiltonian[6,1] = Ω[4] * exp(-im * ϕ[4])

    matrix_hamiltonian[1,8] = Ω[5] * exp(im * ϕ[5])
    matrix_hamiltonian[8,1] = Ω[5] * exp(-im * ϕ[5])

    matrix_hamiltonian[2,5] = Ω[3] * exp(im * ϕ[3])
    matrix_hamiltonian[5,2] = Ω[3] * exp(-im * ϕ[3])

    matrix_hamiltonian[2,7] = Ω[4] * exp(im * ϕ[4])
    matrix_hamiltonian[7,2] = Ω[4] * exp(-im * ϕ[4])

    matrix_hamiltonian[2,9] = Ω[5] * exp(im * ϕ[5])
    matrix_hamiltonian[9,2] = Ω[5] * exp(-im * ϕ[5])

    matrix_hamiltonian[4,5] = Ω[1] * exp(im * ϕ[1])
    matrix_hamiltonian[5,4] = Ω[1] * exp(-im * ϕ[1])

    matrix_hamiltonian[4,10] = Ω[5] * exp(im * ϕ[5])
    matrix_hamiltonian[10,4] = Ω[5] * exp(-im * ϕ[5])

    matrix_hamiltonian[5,11] = Ω[5] * exp(im * ϕ[5])
    matrix_hamiltonian[11,5] = Ω[5] * exp(-im * ϕ[5])

    matrix_hamiltonian[6,7] = Ω[1] * exp(im * ϕ[1])
    matrix_hamiltonian[7,6] = Ω[1] * exp(-im * ϕ[1])

    matrix_hamiltonian[8,9] = Ω[1] * exp(im * ϕ[1])
    matrix_hamiltonian[9,8] = Ω[1] * exp(-im * ϕ[1])

    matrix_hamiltonian[8,10] = Ω[3] * exp(im * ϕ[3])
    matrix_hamiltonian[10,8] = Ω[3] * exp(-im * ϕ[3])

    matrix_hamiltonian[8,10] = Ω[3] * exp(im * ϕ[3])
    matrix_hamiltonian[10,8] = Ω[3] * exp(-im * ϕ[3])

    matrix_hamiltonian[9,11] = Ω[3] * exp(im * ϕ[3])
    matrix_hamiltonian[11,9] = Ω[3] * exp(-im * ϕ[3])

    matrix_hamiltonian[10,11] = Ω[1] * exp(im * ϕ[1])
    matrix_hamiltonian[11,10] = Ω[1] * exp(-im * ϕ[1])

    @test to_matrix(g, Ω, ϕ, Δ) ≈ matrix_hamiltonian
end
