using Lattices
using Test

@testset "Lattices.jl" begin
    for LT in [HoneycombLattice(),
            SquareLattice(), TriangularLattice(),
            LiebLattice(), KagomeLattice(), GeneralLattice(((1.0, 0.0), (0.0, 1.0)), [(0.0, 0.0)])]
        @test generate_sites(LT, 5, 5) |> length == length(latticesites(LT)) * 25
    end
    @test generate_sites(ChainLattice(), 5) |> length == length(latticesites(ChainLattice())) * 5
    l1 = generate_sites(HoneycombLattice(), 5, 5)
    l2 = offsetaxes(l1, -1.0, -2.0)
    l3 = clipaxes(l2, (0.0, 3.0), (0.0, 4.0))
    @test all(loc -> 0 <= loc[1] <= 3 && 0 <= loc[2] <= 4, l3)
    mg = makegrid(l3)
    @test sum(mg.mask) == length(l3) == 14
    @test length(mg.xs) == 6 && length(mg.ys) == 5
    @test length(locations(mg)) == 14
end
