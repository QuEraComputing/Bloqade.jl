using Test, RydbergEmulator
using Random
using Statistics
using LightGraphs, LinearAlgebra

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

@testset "qaoa" begin
    Random.seed!(8)
    g = test_graph()
    params = randn(10)
    ts = params[1:2:end]
    ϕs = params[2:2:end]
    hs = SimpleRydberg.(ϕs)

    # prepair a zero state
    subspace_v = subspace(g)
    st = zeros(ComplexF64, length(subspace_v)); st[1] = 1
    reg = evaluate_qaoa!(SubspaceReg(st, subspace_v), hs, nv(g), ts)
    @test norm(st) ≈ 1
    # 1. sampling
    isets = measure_mis(reg; nshots=10000)
    expected_sampling = Statistics.mean(isets)
    # 2. exact
    expected_exact = expect_mis(reg)
    @test isapprox(expected_exact, expected_sampling; rtol=1e-1)
end
