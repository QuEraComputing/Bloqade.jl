using BloqadeLattices
using Viznet.Compose
using Test, Documenter

@testset "AtomList" begin
    al = AtomList([(0.1, 0.2), (0.3, 0.4), (0.1, 0.8)])
    @test al[2:3] == AtomList([(0.3, 0.4), (0.1, 0.8)])
    @test al[[true, false, true]] == AtomList([(0.1, 0.2), (0.1, 0.8)])
    @test length(al) == 3
    @test size(al) == (3,)
    @test al[3] == (0.1, 0.8)
end

@testset "lattice" begin
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

@testset "neighbors" begin
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
end

@testset "fix site ordering" begin
    lt = generate_sites(KagomeLattice(), 5, 5)
    grd = make_grid(lt[2:end-1])
    x, y = collect_atoms(grd)[1]
    @test issorted(grd.xs)
    @test issorted(grd.ys)
    @test x ≈ 1.0 && y ≈ 0.0
end

@testset "visualize" begin
    lt = generate_sites(KagomeLattice(), 5, 5, scale = 1.5)
    grd = make_grid(lt[2:end-1])
    unitvectors(lattice::AbstractLattice, scale::Real) = [((0.0, 0.0), v .* scale) for v in lattice_vectors(lattice)]
    @test img_atoms(lt; vectors = unitvectors(KagomeLattice(), 1.5)) isa Compose.Context
    # different colors
    @test img_atoms(lt; colors = nothing) isa Compose.Context
    @test img_atoms(lt; colors = "red") isa Compose.Context
    @test img_atoms(lt; colors = fill("blue", length(lt))) isa Compose.Context
    @test img_atoms(lt; colors = ByDensity(randn(length(lt)); vmax = 10)) isa Compose.Context
    @test img_maskedgrid(grd) isa Compose.Context
    @test show(IOBuffer(), MIME"text/html"(), grd) === nothing
    @test show(IOBuffer(), MIME"text/html"(), lt) === nothing
    @test show(IOBuffer(), MIME"image/png"(), grd) === nothing
    @test show(IOBuffer(), MIME"image/png"(), lt) === nothing
end
