using BloqadeLattices
using Test,Distributions,LuxorGraphPlot
using Documenter

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

@testset "cells" begin

    @testset "within_cell" begin
        # 1D Case
        bounds = 4.0
        cell = Parallelepiped(bounds)
        ## Origin
        point = (0.0)
        @test within_cell(cell, point)
        ## Point on line
        point = (1.4)
        @test within_cell(cell, point)
        ## Tip of line
        point = (4.0)
        @test !within_cell(cell, point)
        ## Completely off the line
        point = (11.1)
        @test !within_cell(cell, point)

        # 2D Case
        bounds = zeros((2,2))
        bounds[:,1] .= (3,3)
        bounds[:,2] .= (4,0)
        cell = Parallelepiped(bounds)
        ## origin
        point = (0.0,0.0)
        @test within_cell(cell, point)
        ## lies on one of the accepted sides of the parallelogram
        point = (1.0, 1.0)
        @test within_cell(cell, point)
        ## lies on another accepted side of the parallelogram
        point = (2.0, 0.0)
        @test within_cell(cell, point)
        ## lies inside the parallelogram
        point = (2.0, 1.0)
        @test within_cell(cell, point)
        ## shares a point with an non-accepted side of the parallelogram
        point = (4.0, 0.0)
        @test !within_cell(cell, point)
        ## shares a point with a non-accepted side of the parallelogram
        point = (3.0, 3.0)
        @test !within_cell(cell, point)
        ## lies on the farthest point from the origin and touches two non-accepted sides
        point = (7.0, 3.0)
        @test !within_cell(cell, point)

        # 3D case
        ## Cube
        bounds = zeros((3,3))
        bounds[:,1] .= (1,0,0)
        bounds[:,2] .= (0,1,0)
        bounds[:,3] .= (0,0,1)
        
        cell = Parallelepiped(bounds)
        ## origin
        point = (0.0, 0.0, 0.0)
        @test within_cell(cell, point)
        ## outside parallelpiped
        point = (2.0, 0.0, 0.0)
        @test !within_cell(cell, point)
        ## inside parallelpiped
        point = (0.5, 0.5, 0.5)
        @test within_cell(cell, point)
        ## on a valid surface of the parallelpiped
        point = (0.5, 0.5, 0.0)
        @test within_cell(cell, point)
        ## on a point touching a non-accepted side
        point = (0.0, 0.0, 1.0)
        @test !within_cell(cell, point)
        ## on a non accepted surface
        point = (0.5, 1.0, 0.5)
        @test !within_cell(cell, point)
        ## point touching all non-accepted sides
        point = (1.0, 1.0, 1.0)
        @test !within_cell(cell, point)
    end

    @testset "wrap_around" begin
        # 1D case
        ## Line
        bounds = 3.0
        cell = Parallelepiped(bounds)
        ### points should not wrap around
        points = [0.0, 1.56456, 2.9]
        
        for point in points
            @test all(isapprox.(wrap_around(cell, point),point,atol=1e-15))
        end

        ### points should wrap around
        points = [3.0, 4.0]
        wrapped_points = [0.0, 1.0]

        for (point,wrapped_point) in zip(points,wrapped_points)
            @test all(isapprox.(wrap_around(cell, point),wrapped_point,atol=1e-15))
        end

        # 2D case
        ## Square
        bounds = zeros((2,2))
        bounds[:,1] .= (2,0)
        bounds[:,2] .= (0,2)
        cell = Parallelepiped(bounds)

        ### All of the following points fall outside the square
        ### but should map to its center
        points = [(-1.0, -1.0),
                  (-1.0, 1.0),
                  (-1.0, 3.0),
                  (1.0, 3.0),
                  (3.0, 3.0),
                  (3.0, 1.0),
                  (3.0, -1.0),
                  (1.0, -1.0)]
        for point in points
            @test wrap_around(cell, point) == (1.0, 1.0)
        end

        ### Boundary conditions, points that are on or touch a
        ### side that is not considered part of the cell should be wrapped
        points = [(0.0, 2.0),
                  (2.0, 0.0),
                  (2.0, 2.0)]

        for point in points
            @test wrap_around(cell, point) == (0.0, 0.0)
        end

        ### Randomized testing, generating
        ### points that fall outside bounds of shape, all of which should 
        ### be properly `wrap_around`'ed and
        for (x,y) in zip(rand(Uniform(2.0, 100.0),10), rand(Uniform(2.0, 100.0),10))
            @test within_cell(cell, wrap_around(cell, (x,y)))
        end

        ## Rectangle
        bounds = zeros((2,2))
        bounds[:,1] .= (0.0, 2.0)
        bounds[:,2] .= (4.0, 0.0)
        cell = Parallelepiped(bounds)

        ### All of the following points fall outside the rectangle
        ### but should fall into its interior
        points = [(1.0, -1.0),
                  (1.0, 3.0)]
        
        for point in points
            @test wrap_around(cell, point) == (1.0, 1.0)
        end

        points = [(3.0, -1.0),
                  (3.0, 3.0)]
    
        for point in points
            @test wrap_around(cell, point) == (3.0, 1.0)
        end
        

        ### Boundary conditions, points that are on or touch a
        ### side that is not considered part of the cell should be wrapped
        points = [(0.0, 2.0),
                  (4.0, 0.0),
                  (4.0, 2.0)]
        
        for point in points
            @test wrap_around(cell, point) == (0.0, 0.0)
        end
       
        ### Randomized testing, generating
        ### points that fall outside bounds of shape, all of which should 
        ### be properly `wrap_around`'ed and fall within the shape
        for (x,y) in zip(rand(Uniform(2.0, 100.0),10), rand(Uniform(4.0, 100.0),10))
            @test within_cell(cell, wrap_around(cell, (x,y)))
        end

        ## Parallelogram
        bounds = zeros((2,2))
        bounds[:,1] .= (2,2)
        bounds[:,2] .= (4,0)
        cell = Parallelepiped(bounds)

        ### All of the following points fall outside the parallelogram
        ### but should fall into its interior
        points = [(5.0, 3.0),
                  (-3.5, 0.5),
                  (-5.0, -2.0),
                  (2.0, -1.5)]

        wrapped_points = [(3.0, 1.0),
                          (0.5, 0.5),
                          (1.0, 0.0),
                          (4.0, 0.5)]

        for (point, wrapped_point) in zip(points, wrapped_points)
            @test wrap_around(cell, point) == wrapped_point
        end

        ### Boundary conditions, points that are on or touch a
        ### side that is not considered part of the cell should be wrapped
        points = [(6.0, 2.0),
                  (2.0, 2.0),
                  (4.0, 0.0)]
        
        for point in points
            @test wrap_around(cell, point) == (0.0, 0.0)
        end

        ### Randomized testing, generating
        ### points that fall outside bounds of shape, all of which should 
        ### be properly `wrap_around`'ed and fall within the shape
        for (x,y) in zip(rand(Uniform(6.0, 100.0),10), rand(Uniform(2.0, 100.0),10))
            @test within_cell(cell, wrap_around(cell, (x,y)))
        end

        # 3D Case
        ## Cube
        bounds = zeros((3,3))
        bounds[:,1] .= (2,0,0)
        bounds[:,2] .= (0,2,0)
        bounds[:,3] .= (0,0,2)
        cell = Parallelepiped(bounds)

        ### All of the following points fall outside the parallelogram
        ### but should fall into its interior
        points = [(3.0, 3.0, 3.0)]

        wrapped_points = [(1.0, 1.0, 1.0)]


        for (point, wrapped_point) in zip(points, wrapped_points)
            @test wrap_around(cell, point) == wrapped_point
        end


        ### Boundary conditions, points that are on or touch a
        ### side that is not considered part of the cell should be wrapped
        points = [(2.0, 2.0, 2.0),
                  (0.0, 0.0, 2.0),
                  (0.0, 2.0, 0.0),
                  (2.0, 0.0, 0.0),
                  (2.0, 2.0, 0.0)]
        
        
        for point in points
            @test wrap_around(cell, point) == (0.0, 0.0, 0.0)
        end

        ### Randomized testing, generating
        ### points that fall outside bounds of shape, all of which should 
        ### be properly `wrap_around`'ed and fall within the shape
        for (x,y,z) in zip(rand(Uniform(6.0, 100.0),10), rand(Uniform(2.0, 100.0),10), rand(Uniform(2.0, 100.0), 10))
            @test within_cell(cell, wrap_around(cell, (x,y,z)))
        end
        

    end

    @testset "distance" begin
        # 1D case
        bounds = 3.0
        t = Parallelepiped(bounds)

        point_pairs = [(0.0, 0.1), (1.0, 2.0), (0.5, 2.5), (0.5, 2.1)]
        expected_distances = [0.1, 1.0, 1.0, 1.4]
        for ((x, y),expected_distance) in zip(point_pairs, expected_distances)
            @test isapprox(distance(t, x, y), expected_distance, atol=eps(), rtol=√eps())
        end

        # 2D case
        ## Square
        bounds = zeros((2,2))
        bounds[:,1] .= (3,0)
        bounds[:,2] .= (0,3)
        t = Parallelepiped(bounds)

        point_pairs = [((1.0,1.0), (2.0,2.0)), 
                       ((0.5,0.5), (2.5,2.5)),
                       ((0.4,2.4), (2.0,2.4)),
                       ((2.5,2.5), (2.0,0.5))]
        expected_distances = [√2.0, 2√0.5, 1.4, √1.25]
        for ((x, y),expected_distance) in zip(point_pairs, expected_distances)
            @test isapprox(distance(t, x, y), expected_distance, atol=eps(), rtol=√eps())
        end

        ## Rectangle
        bounds = zeros((2,2))
        bounds[:,1] .= (0,2)
        bounds[:,2] .= (4,0)
        t = Parallelepiped(bounds)

        point_pairs = [((1.0, 1.0), (2.0, 1.0)),
                       ((3.0, 0.4), (3.0, 1.5)),
                       ((2.0, 0.0), (3.5, 1.5))]
        expected_distances = [1.0, 0.9, √2.5]
        for ((x, y),expected_distance) in zip(point_pairs, expected_distances)
            @test isapprox(distance(t, x, y), expected_distance, atol=eps(), rtol=√eps())
        end

        ## Parallelogram
        bounds = zeros((2,2))
        bounds[:,1] .= (3,3)
        bounds[:,2] .= (4,0)
        t = Parallelepiped(bounds)

        point_pairs = [((4.0, 1.0), (5.0, 2.0)),
                       ((1.0, 0.5), (6.0, 2.5)),
                       ((3.0, 1.0), (5.0, 2.0))]
        expected_distances = [√2.0, 2√1.25, √5.0]
        for ((x, y),expected_distance) in zip(point_pairs, expected_distances)
            @test isapprox(distance(t, x, y), expected_distance, atol=eps(), rtol=√eps())
        end

        # 3D
        ## Cube
        bounds = zeros((3,3))
        bounds[:,1] .= (3,0,0)
        bounds[:,2] .= (0,3,0)
        bounds[:,3] .= (0,0,3)
        t = Parallelepiped(bounds)

        point_pairs = [((1.0, 1.0, 1.0), (2.0, 2.0, 2.0)),
                       ((2.5, 2.5, 2.9), (2.5, 2.5, 1.1))]
        expected_distances = [√3, 1.2]
        for ((x, y),expected_distance) in zip(point_pairs, expected_distances)
            @test isapprox(distance(t, x, y), expected_distance, atol=eps(), rtol=√eps())
        end

    end

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
    BloqadeLattices.darktheme!()
    lt = generate_sites(KagomeLattice(), 5, 5, scale = 1.5)
    grd = make_grid(lt[2:end-1])
    unitvectors(lattice::AbstractLattice, scale::Real) = [((0.0, 0.0), v .* scale) for v in lattice_vectors(lattice)]
    @test img_atoms(lt; vectors = unitvectors(KagomeLattice(), 1.5)) isa LuxorGraphPlot.Drawing
    # different colors
    @test img_atoms(lt; colors = nothing) isa LuxorGraphPlot.Drawing
    @test img_atoms(lt; node_fill_color = "red") isa LuxorGraphPlot.Drawing
    @test img_atoms(lt; colors = fill("blue", length(lt))) isa LuxorGraphPlot.Drawing
    @test img_atoms(lt; colors = ByDensity(randn(length(lt)); vmax = 10)) isa LuxorGraphPlot.Drawing
    @test img_maskedgrid(grd) isa LuxorGraphPlot.Drawing
    @test show(IOBuffer(), MIME"image/svg+xml"(), grd) === nothing
    @test show(IOBuffer(), MIME"image/svg+xml"(), lt) === nothing
    @test show(IOBuffer(), MIME"image/png"(), grd) === nothing
    @test show(IOBuffer(), MIME"image/png"(), lt) === nothing

    BloqadeLattices.lighttheme!()
    @test img_atoms(lt; colors = nothing) isa LuxorGraphPlot.Drawing
    @test show(IOBuffer(), MIME"image/svg+xml"(), lt) === nothing
end