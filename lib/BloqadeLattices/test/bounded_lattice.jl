using Test
using BloqadeLattices


for scale in [1,1.0,2,2.0,π,1//2]
    # chain
    bounded_lattice = parallelepiped_region(ChainLattice(),(4,);pbc=true,scale=scale)
    expected_positions = [(0.0*scale,),(1.0*scale,),(2.0*scale,),(3.0*scale,)]
    @test issetequal(bounded_lattice.site_positions,expected_positions)
    @test dimension(bounded_lattice) == 1

    # square
    bounded_lattice = parallelepiped_region(SquareLattice(),(2,0),(0,2);pbc=true,scale=scale)
    expected_positions = [(0.0*scale,0.0*scale),(1.0*scale,0.0*scale),(0.0*scale,1.0*scale),(1.0*scale,1.0*scale)]
    @test issetequal(bounded_lattice.site_positions,expected_positions)
    @test dimension(bounded_lattice) == 2

    # tilted square
    bounded_lattice = parallelepiped_region(SquareLattice(),(3,2),(-2,3);scale=scale)
    expected_positions = [
        (-1.0*scale, 2.0*scale), (-1.0*scale, 3.0*scale), (0.0*scale, 0.0*scale), (0.0*scale, 1.0*scale), 
        (0.0*scale, 2.0*scale), (0.0*scale, 3.0*scale), (0.0*scale, 4.0*scale), (1.0*scale, 1.0*scale), 
        (1.0*scale, 2.0*scale), (1.0*scale, 3.0*scale), (1.0*scale, 4.0*scale), (2.0*scale, 2.0*scale), (2.0*scale, 3.0*scale)
    ]
    @test issetequal(bounded_lattice.site_positions,expected_positions)

    # kagome rectangle boundary
    bounded_lattice = parallelepiped_region(KagomeLattice(),(2,2),(-2,2))
    expected_positions = [
        (-0.75*scale, 1.299038105676658*scale), (-0.5*scale, 0.8660254037844386*scale), (-0.25*scale, 0.4330127018922193*scale), 
        (-0.25*scale, 1.299038105676658*scale), (0.0*scale, 0.0*scale), (0.0*scale, 1.7320508075688772*scale), 
        (0.25*scale, 0.4330127018922193*scale), (0.25*scale, 1.299038105676658*scale), (0.25*scale, 2.1650635094610964*scale), 
        (0.5*scale, 0.8660254037844386*scale), (0.75*scale, 0.4330127018922193*scale), (0.75*scale, 1.299038105676658*scale), 
        (0.75*scale, 2.1650635094610964*scale), (1.0*scale, 1.7320508075688772*scale), (1.25*scale, 1.299038105676658*scale), 
        (1.25*scale, 2.1650635094610964*scale), (1.5*scale, 0.8660254037844386*scale), (1.5*scale, 2.598076211353316*scale), 
        (1.75*scale, 1.299038105676658*scale), (1.75*scale, 2.1650635094610964*scale), (1.75*scale, 3.031088913245535*scale), 
        (2.0*scale, 1.7320508075688772*scale), (2.25*scale, 1.299038105676658*scale), (2.25*scale, 2.1650635094610964*scale)
    ]
    @test issetequal(bounded_lattice.site_positions,expected_positions)
    @test dimension(bounded_lattice) == 2

    # cube
    cube_vectors = ((1.0, 0.0, 0.0), 
                    (0.0, 1.0, 0.0), 
                    (0.0, 0.0, 1.0))
    cube_sites = ((0.0, 0.0, 0.0),)
    cube_lattice = GeneralLattice(cube_vectors, cube_sites)
    bounded_lattice = parallelepiped_region(cube_lattice, 
                                            (2, 0, 0), 
                                            (0, 2, 0),
                                            (0, 0, 2);
                                            pbc=false,scale=scale)
    expected_positions = [
        (0.0*scale, 0.0*scale, 0.0*scale), 
        (1.0*scale, 0.0*scale, 0.0*scale), (0.0, 1.0*scale, 0.0*scale), (0.0*scale, 0.0*scale, 1.0*scale),
        (1.0*scale, 1.0*scale, 0.0*scale), (1.0*scale, 0.0*scale, 1.0*scale), (0.0*scale, 1.0*scale, 1.0*scale),
        (1.0*scale, 1.0*scale, 1.0*scale)
    ]
    @test issetequal(bounded_lattice.site_positions, expected_positions)
end