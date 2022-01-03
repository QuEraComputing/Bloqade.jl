using Test
using Random
using EaRydCore
using Graphs

@testset "unit_disk" begin
    g = SimpleGraph(4)
    positions = [0 0 0.3 0.4; 0 0.6 0 0.6]
    atms = rydatoms(positions)
    g2 = unit_disk_graph(atms)
    @test ne(g2) == 6

    radius = 0.5
    for edge = [(1,3),(2,4)]
        add_edge!(g,edge)
    end

    h = unit_disk_graph(atms,radius)
    @test g == h

    radius = sqrt(0.3^2 + 0.6^2) + 0.01
    g2 = complete_graph(4)
    rem_edge!(g2,1,4)

    h2 = unit_disk_graph(atms,radius)
    @test g2 == h2

    Random.seed!(42)
    g1 = unit_disk_graph(rand_atoms(5, 0.8; ndims=2), 1.0)

    Random.seed!(42)
    g2 = rand_unit_disk_graph(5, 0.8)

    @test g1 == g2
end
