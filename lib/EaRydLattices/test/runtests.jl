using EaRydLattices
using Viznet.Compose
using Test

@testset "lattice" begin
    for LT in [HoneycombLattice(),
            SquareLattice(), TriangularLattice(),
            LiebLattice(), KagomeLattice(), GeneralLattice(((1.0, 0.0), (0.0, 1.0)), [(0.0, 0.0)])]
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
    
    # rescale axes
    sites = AtomList([(0.2, 0.3), (0.4, 0.8)])
    @test (sites |> rescale_axes(2.0)) == [(0.4, 0.6), (0.8, 1.6)]
end

@testset "neighbors" begin
    for (lt, i, s1, s2) in [(HoneycombLattice(), 25, 3, 6),
            (SquareLattice(), 13, 4, 4), (TriangularLattice(), 13, 6, 6),
            (LiebLattice(), 37, 4, 4), (KagomeLattice(), 37, 4, 4), (GeneralLattice(((1.0, 0.0), (0.0, 1.0)), [(0.0, 0.0)]), 13, 4, 4)]
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
    lt = generate_sites(KagomeLattice(), 5, 5)
    grd = make_grid(lt[2:end-1])
    @test viz_atoms(IOBuffer(), lt) === nothing
    @test viz_maskedgrid(IOBuffer(), grd) === nothing
    @test show(IOBuffer(), MIME"text/html"(), KagomeLattice()) === nothing
    @test show(IOBuffer(), MIME"text/html"(), ChainLattice()) === nothing
    @test show(IOBuffer(), MIME"text/html"(), grd) === nothing
    @test show(IOBuffer(), MIME"text/html"(), lt) === nothing
end