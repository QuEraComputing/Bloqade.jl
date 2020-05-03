using Test, RydbergEmulator
using Random
using LightGraphs
using Evolutionary

"""
    cmaes_train_mis(graph, ϕs0, ts0)

Obtain the MIS using CMA-ES training. Return the optimal `ϕs` and `ts`.
"""
function cmaes_train_mis(graph, ϕs0, ts0)
    @assert length(ϕs0) == length(ts0)
    p = length(ϕs0)
    params = vcat(ϕs0, ts0)
    res = Evolutionary.optimize(params, CMAES(μ = 5, λ = 100)) do params
        p = length(params)÷2
        ϕs = params[1:p]
        ts = params[p+1:end]
        -mean_independent_set(graph, ϕs, ts; nshots=nothing)
    end
    res.minimizer[1:p], res.minimizer[p+1:end]
end

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
