using Lattices
using Test

@testset "Lattices.jl" begin
    l1 = generate_sites(HoneycombLattice(), 5, 5)
    l2 = offset(l1, -1.0, -2.0)
    l3 = clip(l2, (-1.0, 3.0), (-2.0, 4.0))
    @test all(loc -> -1 <= loc[1] <= 3 && -2 <= loc[2] <= 4, l3)
end
