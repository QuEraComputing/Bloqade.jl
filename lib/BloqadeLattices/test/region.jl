using Test
using BloqadeLattices

@testset "within_region" begin
    # 1D Case
    bounds = 4.0
    region = Parallelepiped(bounds)
    ## Origin
    point = (0.0)
    @test point ∈ region
    ## Point on line
    point = (1.4)
    @test point ∈ region
    ## Tip of line
    point = (4.0)
    @test !(point ∈ region)
    ## Completely off the line
    point = (11.1)
    @test !(point ∈ region)

    # 2D Case
    bounds = zeros((2,2))
    bounds[:,1] .= (3,3)
    bounds[:,2] .= (4,0)
    region = Parallelepiped(bounds)
    ## origin
    point = (0.0,0.0)
    @test point ∈ region
    ## lies on one of the accepted sides of the parallelogram
    point = (1.0, 1.0)
    @test point ∈ region
    ## lies on another accepted side of the parallelogram
    point = (2.0, 0.0)
    @test point ∈ region
    ## lies inside the parallelogram
    point = (2.0, 1.0)
    @test point ∈ region
    ## shares a point with an non-accepted side of the parallelogram
    point = (4.0, 0.0)
    @test !(point ∈ region)
    ## shares a point with a non-accepted side of the parallelogram
    point = (3.0, 3.0)
    @test !(point ∈ region)
    ## lies on the farthest point from the origin and touches two non-accepted sides
    point = (7.0, 3.0)
    @test !(point ∈ region)

    # 3D case
    ## Cube
    bounds = zeros((3,3))
    bounds[:,1] .= (1,0,0)
    bounds[:,2] .= (0,1,0)
    bounds[:,3] .= (0,0,1)
    
    region = Parallelepiped(bounds)
    ## origin
    point = (0.0, 0.0, 0.0)
    @test point ∈ region
    ## outside parallelpiped
    point = (2.0, 0.0, 0.0)
    @test !(point ∈ region)
    ## inside parallelpiped
    point = (0.5, 0.5, 0.5)
    @test point ∈ region
    ## on a valid surface of the parallelpiped
    point = (0.5, 0.5, 0.0)
    @test point ∈ region
    ## on a point touching a non-accepted side
    point = (0.0, 0.0, 1.0)
    @test !(point ∈ region)
    ## on a non accepted surface
    point = (0.5, 1.0, 0.5)
    @test !(point ∈ region)
    ## point touching all non-accepted sides
    point = (1.0, 1.0, 1.0)
    @test !(point ∈ region)
end

@testset "mod" begin

    ## Use this to generate random values from uniform distribution
    uniform(a,b) = (a-b)*rand() + b

    # 1D case
    ## Line
    bounds = 3.0
    region = Parallelepiped(bounds)
    ### points should not wrap around
    points = [0.0, 1.56456, 2.9]
    
    for point in points
        @test all(isapprox.(mod(point,region),point,atol=1e-15))
    end

    ### points should wrap around
    points = [3.0, 4.0]
    wrapped_points = [0.0, 1.0]

    for (point,wrapped_point) in zip(points,wrapped_points)
        @test all(isapprox.(mod(point,region),wrapped_point,atol=1e-15))
    end

    ### Randomized testing, generating points
    ### that fall outside the line which should,
    ### after wrapping, fall inside the line
    for x in [uniform(-100.0, 100.0) for _ in 1:10]
        @test mod(x,region) ∈ region
    end

    # 2D case
    ## Square
    bounds = zeros((2,2))
    bounds[:,1] .= (2,0)
    bounds[:,2] .= (0,2)
    region = Parallelepiped(bounds)

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
        @test mod(point,region) == (1.0, 1.0)
    end

    ### Boundary conditions, points that are on or touch a
    ### side that is not considered part of the region should be wrapped
    points = [(0.0, 2.0),
              (2.0, 0.0),
              (2.0, 2.0)]

    for point in points
        @test mod(point,region) == (0.0, 0.0)
    end


    ### Randomized testing, generating points
    ### that fall outside the square which should,
    ### after wrapping, fall inside the square
    for (x,y) in zip([uniform(-100.0, 100.0) for _ in 1:10], 
                     [uniform(-100.0, 100.0) for _ in 1:10])
        @test mod((x,y),region) ∈ region
    end

    ## Rectangle
    bounds = zeros((2,2))
    bounds[:,1] .= (0.0, 2.0)
    bounds[:,2] .= (4.0, 0.0)
    region = Parallelepiped(bounds)

    ### All of the following points fall outside the rectangle
    ### but should fall into its interior
    points = [(1.0, -1.0),
              (1.0, 3.0)]
    
    for point in points
        @test mod(point,region) == (1.0, 1.0)
    end

    points = [(3.0, -1.0),
              (3.0, 3.0)]

    for point in points
        @test mod(point,region) == (3.0, 1.0)
    end
    

    ### Boundary conditions, points that are on or touch a
    ### side that is not considered part of the region should be wrapped
    points = [(0.0, 2.0),
              (4.0, 0.0),
              (4.0, 2.0)]
    
    for point in points
        @test mod(point,region) == (0.0, 0.0)
    end
   

    ### Randomized testing, generating points
    ### that fall outside the rectangle which should,
    ### after wrapping, fall inside the rectangle
    for (x,y) in zip([uniform(-100.0, 100.0) for _ in 1:10], 
                     [uniform(-100.0, 100.0) for _ in 1:10])
        @test mod((x,y),region) ∈ region
    end

    ## Parallelogram
    bounds = zeros((2,2))
    bounds[:,1] .= (2,2)
    bounds[:,2] .= (4,0)
    region = Parallelepiped(bounds)

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
        @test mod(point,region) == wrapped_point
    end

    ### Boundary conditions, points that are on or touch a
    ### side that is not considered part of the region should be wrapped
    points = [(6.0, 2.0),
              (2.0, 2.0),
              (4.0, 0.0)]
    
    for point in points
        @test mod(point,region) == (0.0, 0.0)
    end

    ### Randomized testing, generating points
    ### that fall outside the parallelogram which should,
    ### after wrapping, fall inside the square
    for (x,y) in zip([uniform(-100.0, 100.0) for _ in 1:10], 
                     [uniform(-100.0, 100.0) for _ in 1:10])
        @test mod((x,y),region) ∈ region
    end

    # 3D Case
    ## Cube
    bounds = zeros((3,3))
    bounds[:,1] .= (2,0,0)
    bounds[:,2] .= (0,2,0)
    bounds[:,3] .= (0,0,2)
    region = Parallelepiped(bounds)

    ### All of the following points fall outside the parallelogram
    ### but should fall into its interior
    points = [(3.0, 3.0, 3.0)]

    wrapped_points = [(1.0, 1.0, 1.0)]


    for (point, wrapped_point) in zip(points, wrapped_points)
        @test mod(point,region) == wrapped_point
    end


    ### Boundary conditions, points that are on or touch a
    ### side that is not considered part of the region should be wrapped
    points = [(2.0, 2.0, 2.0),
              (0.0, 0.0, 2.0),
              (0.0, 2.0, 0.0),
              (2.0, 0.0, 0.0),
              (2.0, 2.0, 0.0)]
    
    
    for point in points
        @test mod(point,region) == (0.0, 0.0, 0.0)
    end

    ### Randomized testing, generating points
    ### that fall outside the cube which should,
    ### after wrapping, fall inside the cube
    for (x,y,z) in zip([uniform(-100.0, 100.0) for _ in 1:10], 
                       [uniform(-100.0, 100.0) for _ in 1:10], 
                       [uniform(-100.0, 100.0) for _ in 1:10])
        @test mod((x,y,z),region) ∈ region
    end
end

@testset "distance" begin
    # 1D case
    bounds = 3.0
    t = Parallelepiped(bounds)

    point_pairs = [(0.0, 0.1), (1.0, 2.0), (0.5, 2.5), (0.5, 2.1),(2.9,0.1)]
    expected_distances = [0.1, 1.0, 1.0, 1.4, 0.2]
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
    expected_distances = [√2.0, √2, 1.4, √1.25]
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
    expected_distances = [√2.0, 2*√1.25, √5.0]
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

@testset "generate_sites_in_region" begin

    # Chain Lattice
    lattice = ChainLattice()
    region = Parallelepiped(4.0)
    expected_positions = [(0.0,), (1.0,), (2.0,), (3.0,)]
    @test issetequal(generate_sites_in_region(lattice, region), expected_positions)
    ## Negative bounds
    region = Parallelepiped(-4.0)
    expected_positions = [(0.0,), (-1.0,), (-2.0,), (-3.0,)]
    @test issetequal(generate_sites_in_region(lattice, region), expected_positions)

    # Square Lattice
    lattice = SquareLattice()
    bounds = zeros((2,2))
    bounds[:,1] .= (0.0, 2.0)
    bounds[:,2] .= (2.0, 0.0)
    region = Parallelepiped(bounds)
    expected_positions = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0)]
    @test issetequal(generate_sites_in_region(lattice,region), expected_positions)
    ## Negative bounds
    bounds[:,1] .= (2.0, 2.0)
    bounds[:,2] .= (2.0, -2.0)
    expected_positions = [(0,0), (1,1), (1,0), (1,-1), (2,1), (2,0), (2,-1), (3,0)]
    region = Parallelepiped(bounds)
    @test issetequal(generate_sites_in_region(lattice, region), expected_positions)
    
    # Triangular Lattice
    lattice = TriangularLattice()
    bounds = zeros((2,2))
    bounds[:,1] .= (0, 2)
    bounds[:,2] .= (2, 0)
    region = Parallelepiped(bounds)
    expected_positions = [(0.0, 1.7320508075688772), 
                      (0.0, 0.0), 
                      (0.5, 0.8660254037844386), 
                      (1.0, 1.7320508075688772), 
                      (1.0, 0.0), 
                      (1.5, 0.8660254037844386)]
    @test issetequal(generate_sites_in_region(lattice, region), expected_positions)

    # generate rotation matrix
    rot_mat(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]

    # Honeycomb Lattice
    lattice = HoneycombLattice()
    bounds = zeros((2,2))
    # create region that slices through lattice unit cells.
    # Start off with vectors for Honeycomb lattice, 
    # then rotate and scale as needed
    (a1,a2) = lattice_vectors(lattice)      
    bounds[:,1] .= rot_mat(deg2rad(50)) * [a1...,] * 2
    bounds[:,2] .= rot_mat(deg2rad(50)) * [a2...,] * 2

    region = Parallelepiped(bounds)

    expected_positions = [
        (-0.5, 2.0207259421636903), 
        (0.0, 0.0), 
        (0.0, 1.1547005383792515), 
        (0.0, 1.7320508075688772), 
        (0.5, 0.8660254037844386), 
        (0.5, 2.0207259421636903), 
        (0.5, 2.598076211353316), 
        (1.0, 1.7320508075688772)
    ]
    @test issetequal(generate_sites_in_region(lattice, region), expected_positions)

end