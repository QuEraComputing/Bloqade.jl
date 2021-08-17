using LightGraphs
using Yao
using ContinuousEmulator
using Test

test_graph = SimpleGraph(5)
add_edge!(test_graph, 1, 2)
add_edge!(test_graph, 2, 3)
add_edge!(test_graph, 2, 4)
add_edge!(test_graph, 2, 5)
add_edge!(test_graph, 3, 4)
add_edge!(test_graph, 4, 5)

test_subspace_v = [0, 1, 2, 4, 5, 8, 9, 16, 17, 20, 21]
test_subspace = blockade_subspace(test_graph)

@testset "contiguous time" begin
    h = XTerm(5, 1.0, sin)
    dt = 1e-5
    @testset "subspace" begin
        r1 = RydbergEmulator.zero_state(5, test_subspace)
        emulate!(r1, 0.2, h)

        r2 = RydbergEmulator.zero_state(5, test_subspace)
        emulate!(r2, map(_->dt, 0.0:dt:0.2), map(h, 0.0:dt:0.2))

        @test isapprox(r1.state, r2.state; atol=1e-4)
    end

    @testset "fullspace" begin
        r1 = Yao.zero_state(5)
        emulate!(r1, 0.2, h)
        r2 = Yao.zero_state(5)
        emulate!(r2, map(_->dt, 0.0:dt:0.2), map(h, 0.0:dt:0.2))
        @test isapprox(r1.state, r2.state; atol=1e-4)
    end
end
