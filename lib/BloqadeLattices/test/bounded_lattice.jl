using Test
using BloqadeLattices

# chain
bounded_lattice = parallelepiped_region(ChainLattice(),(4,);pbc=true)
expected_positions = [(0.0,),(1.0,),(2.0,),(3.0,)]
@test issetequal(bounded_lattice.site_positions,expected_positions)
@test dimension(bounded_lattice) == 1

# square
bounded_lattice = parallelepiped_region(SquareLattice(),(2,0),(0,2);pbc=true)
expected_positions = [(0.0,0.0),(1.0,0.0),(0.0,1.0),(1.0,1.0)]
@test issetequal(bounded_lattice.site_positions,expected_positions)
@test dimension(bounded_lattice) == 2

# tilted square
bounded_lattice = parallelepiped_region(SquareLattice(),(3,2),(-2,3))
expected_positions = [
    (-1.0, 2.0), (-1.0, 3.0), (0.0, 0.0), (0.0, 1.0), 
    (0.0, 2.0), (0.0, 3.0), (0.0, 4.0), (1.0, 1.0), 
    (1.0, 2.0), (1.0, 3.0), (1.0, 4.0), (2.0, 2.0), (2.0, 3.0)
]
@test issetequal(bounded_lattice.site_positions,expected_positions)

# kagome rectangle boundary
bounded_lattice = parallelepiped_region(KagomeLattice(),(2,2),(-2,2))
expected_positions = [
    (-0.75, 1.299038105676658), (-0.5, 0.8660254037844386), (-0.25, 0.4330127018922193), 
    (-0.25, 1.299038105676658), (0.0, 0.0), (0.0, 1.7320508075688772), 
    (0.25, 0.4330127018922193), (0.25, 1.299038105676658), (0.25, 2.1650635094610964), 
    (0.5, 0.8660254037844386), (0.75, 0.4330127018922193), (0.75, 1.299038105676658), 
    (0.75, 2.1650635094610964), (1.0, 1.7320508075688772), (1.25, 1.299038105676658), 
    (1.25, 2.1650635094610964), (1.5, 0.8660254037844386), (1.5, 2.598076211353316), 
    (1.75, 1.299038105676658), (1.75, 2.1650635094610964), (1.75, 3.031088913245535), 
    (2.0, 1.7320508075688772), (2.25, 1.299038105676658), (2.25, 2.1650635094610964)
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
                                        pbc=false)
expected_positions = [
    (0.0, 0.0, 0.0), 
    (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0),
    (1.0, 1.0, 0.0), (1.0, 0.0, 1.0), (0.0, 1.0, 1.0),
    (1.0, 1.0, 1.0)
]
@test issetequal(bounded_lattice.site_positions, expected_positions)