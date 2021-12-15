# Copyright 2021 QuEra Computing Inc. All rights reserved.

module EaRydLattices

using NearestNeighbors
using Viznet: Viznet
using Viznet.Compose

export # types
    AbstractLattice, BravaisLattice, HoneycombLattice,
    SquareLattice, TriangularLattice, ChainLattice,
    LiebLattice, KagomeLattice, GeneralLattice,
    # interfaces
    bravais, generate_sites, offset_axes,
    clip_axes, lattice_sites, lattice_vectors,
    # grid
    MaskedGrid, make_grid, locations, random_dropout,
    # nearest neighbors
    make_kdtree, grouped_nearest,
    # visualize
    viz_atoms, viz_maskedgrid

include("lattice.jl")
include("neighbors.jl")
include("visualize.jl")

end
