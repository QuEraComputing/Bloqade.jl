using Test
using RydbergEmulator
using LightGraphs
using OrderedCollections

include("utils.jl")

@testset "simple graph hamiltonian subspace" begin
    g = SimpleGraph(5)
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)
    add_edge!(g, 2, 4)
    add_edge!(g, 2, 5)
    add_edge!(g, 3, 4)
    add_edge!(g, 4, 5)

    subspace = Subspace(g)
    set_subspace = [0,1,2,4,8,16,5,9,17,20,21]
    @test collect(keys(subspace)) == sort(set_subspace)
    @test collect(values(subspace)) == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

    # test to_matrix (hamiltonian creation)
    Ω = Float64[2.5, 3.4, 0.2, 1.7, 4.3]
    Δ = Float64[1.2, 3.4, 2.6, 0.2, 1.8]
    ϕ = Float64[0.5, 0.4, -0.2, -1.2, 10.2]

    target = create_test_hamiltonian(Δ, Ω, ϕ)
    h = XTerm(Ω, ϕ) + ZTerm(Δ)
    @test to_matrix(h, subspace) ≈ target

    H = to_matrix(h, subspace)
    @test H ≈ update_term!(copy(H), h, subspace)
end
