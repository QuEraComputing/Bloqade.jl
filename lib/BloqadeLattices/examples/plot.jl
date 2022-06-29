using BloqadeLattices

function ring12()
    nsites = 12; # 12-site ring
    unit_disk_radius = 6.9 # Distance between nearest-neighbor atoms, in microns; R_min

    R = unit_disk_radius / (2 * sin(2 * pi / (nsites) / 2)) # Radius of the circle, using a little trigonometry; it is also the next-nearest neighbor distance, R_max.
    pos = [(R * sin(i * 2 * pi / (nsites)), R * cos(i * 2 * pi / (nsites))) for i in 1:nsites] # Positions of each atom
    AtomList(pos); # Define the atom positions as an AtomList.
end

unitvectors(lattice::AbstractLattice{2}) = [((0.0, 0.0), v) for v in lattice_vectors(lattice)]
BloqadeLattices.img_atoms(ring12())
square = SquareLattice()
BloqadeLattices.img_atoms(generate_sites(square, 10, 10); vectors=unitvectors(square), bond_linewidth=0.015)
honeycomb = HoneycombLattice()
BloqadeLattices.img_atoms(generate_sites(honeycomb, 5, 5); vectors=unitvectors(honeycomb), bond_linewidth=0.015)