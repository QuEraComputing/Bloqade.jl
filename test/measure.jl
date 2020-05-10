using Test, RydbergEmulator
using Random
using LightGraphs, LinearAlgebra
using BitBasis

function test_graph()
    g = SimpleGraph(5)
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)
    add_edge!(g, 2, 4)
    add_edge!(g, 2, 5)
    add_edge!(g, 3, 4)
    add_edge!(g, 4, 5)
    return g
end

@testset "test measurement" begin
    Random.seed!(8)
    g = test_graph()
    params = randn(10)
    ts = params[1:2:end]
    ϕs = params[2:2:end]
    hs = SimpleRydberg.(ϕs)

    # prepare a zero state
    subspace_v = subspace(g)
    r = RydbergEmulator.zero_state(5, subspace_v)
    qaoa = QAOA{5}(subspace_v, hs, ts)
    r |> qaoa
    Random.seed!(5)
    # 1. sampling
    samples = measure(r; nshots=10000)
    @test samples isa Vector{<:BitStr}
    expected_sampling = samples .|> count_vertices |> mean
    # 2. exact
    expected_exact = r |> mean_nv
    @test isapprox(expected_exact, expected_sampling; rtol=1e-1)

    ####### softmis
    Random.seed!(5)
    # 1. sampling
    expected_sampling = r |> measure(nshots=10000) |> soft_misloss(0.5)
    # 2. exact
    expected_exact = r |> soft_misloss(0.5)
    @test isapprox(expected_exact, expected_sampling; rtol=1e-1)
end
