using Test
using RydbergEmulator
using LightGraphs

@testset "unit_disk" begin
    g = SimpleGraph(4)
    positions = [0 0 0.3 0.4; 0 0.6 0 0.6]
    atms = rydatoms(positions)

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
end
