using Test
using BloqadeLattices

for (lt, i, s1, s2) in [
    (HoneycombLattice(), 25, 3, 6),
    (SquareLattice(), 13, 4, 4),
    (TriangularLattice(), 13, 6, 6),
    (LiebLattice(), 37, 4, 4),
    (KagomeLattice(), 37, 4, 4),
    (GeneralLattice(((1.0, 0.0), (0.0, 1.0)), [(0.0, 0.0)]), 13, 4, 4),
]
    atoms = generate_sites(lt, 5, 5)
    tree = make_kdtree(atoms)
    @show lt
    gn = grouped_nearest(tree, i, 20)
    @test gn[0] == [i]
    @test length(gn[1]) == s1
    @test length(gn[2]) == s2
    @test length(gn[10]) == 0
end