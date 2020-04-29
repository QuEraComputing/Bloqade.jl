using Test, RydbergEmulator
using Random
using LightGraphs

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

@testset "mis" begin
    Random.seed!(2)
    graph = test_graph()
    p = 5
    @test isapprox(mean_independent_set(graph, rand(p), zeros(p); nshots=nothing), 0, atol=1e-8)
    ϕs0 = rand(p) * 2π
    ts0 = rand(p)
    ϕs, ts = cmaes_train_mis(graph, ϕs0, ts0)
    n = mean_independent_set(graph, ϕs, ts; nshots=nothing)
    @test n > mean_independent_set(graph, ϕs0, ts0)
    @test isapprox(n, 3, atol=0.25)

    reg = qaoa_on_graph(graph, ϕs::AbstractVector, ts::AbstractVector)
    res = measure(reg, nshots=100)
    @test count(x->(x==21), res) > 50
end
