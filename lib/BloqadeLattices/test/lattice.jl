using BloqadeLattices
using Test

@testset "AtomList" begin
    al = AtomList([(0.1, 0.2), (0.3, 0.4), (0.1, 0.8)])
    @test al[2:3] == AtomList([(0.3, 0.4), (0.1, 0.8)])
    @test al[[true, false, true]] == AtomList([(0.1, 0.2), (0.1, 0.8)])
    @test length(al) == 3
    @test size(al) == (3,)
    @test al[3] == (0.1, 0.8)
end

@testset "lattices" begin
    for LT in [
        HoneycombLattice(),
        SquareLattice(),
        TriangularLattice(),
        LiebLattice(),
        KagomeLattice(),
        GeneralLattice(((1.0, 0.0), (0.0, 1.0)), [(0.0, 0.0)]),
    ]
        @test BloqadeLattices.dimension(LT) == 2
        @test generate_sites(LT, 5, 5) |> length == length(lattice_sites(LT)) * 25
    end
    lt1 = generate_sites(ChainLattice(), 5)
    @test lt1 |> length == length(lattice_sites(ChainLattice())) * 5
    @test make_grid(lt1) isa MaskedGrid
    l1 = generate_sites(HoneycombLattice(), 5, 5)
    l2 = l1 |> offset_axes(-1.0, -2.0)
    l3 = l2 |> clip_axes((0.0, 3.0), (0.0, 4.0))
    @test all(loc -> 0 <= loc[1] <= 3 && 0 <= loc[2] <= 4, l3)
    mg = make_grid(l3)
    @test sum(mg.mask) == length(l3) == 14
    @test length(mg.xs) == 6 && length(mg.ys) == 5
    @test length(collect_atoms(mg)) == 14
    l4 = l3 |> random_dropout(1.0)
    @test length(l4) == 0
    l4 = random_dropout(l3, 0.0)
    @test length(l4) == length(l3)
    l5 = random_dropout(l3, 0.5)
    @test length(l5) == 7
    @test_throws ArgumentError random_dropout(l3, -0.5)

    # Rectangular Lattice Defaults
    rectangular_lattice = RectangularLattice(1.0)
    @test lattice_sites(rectangular_lattice) == ((0.0, 0.0),)
    @test lattice_vectors(rectangular_lattice)[2][2] == rectangular_lattice.aspect_ratio

    # Lattice Dimensions

    # rescale axes
    sites = AtomList([(0.2, 0.3), (0.4, 0.8)])
    @test (sites |> rescale_axes(2.0)) == [(0.4, 0.6), (0.8, 1.6)]

end

@testset "fix site ordering" begin
    lt = generate_sites(KagomeLattice(), 5, 5)
    grd = make_grid(lt[2:end-1])
    x, y = collect_atoms(grd)[1]
    @test issorted(grd.xs)
    @test issorted(grd.ys)
    @test x ≈ 1.0 && y ≈ 0.0
end