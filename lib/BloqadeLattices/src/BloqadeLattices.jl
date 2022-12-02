# Copyright 2021 QuEra Computing Inc. All rights reserved.

module BloqadeLattices

using NearestNeighbors
using StatsBase
using LuxorGraphPlot
using LuxorGraphPlot: Point
using LuxorGraphPlot.Luxor: Colors
using LinearAlgebra

export # types
    AbstractLattice,
    HoneycombLattice,
    SquareLattice,
    TriangularLattice,
    ChainLattice,
    LiebLattice,
    KagomeLattice,
    GeneralLattice,
    RectangularLattice,
    # interfaces
    generate_sites,
    offset_axes,
    random_dropout,
    rescale_axes,
    clip_axes,
    lattice_sites,
    lattice_vectors,
    # grid
    MaskedGrid,
    make_grid,
    AtomList,
    collect_atoms,
    # nearest neighbors
    make_kdtree,
    grouped_nearest,
    DistanceGroup,
    # visualize
    img_atoms,
    img_maskedgrid,
    ByDensity,
    # bounds
    AbstractRegion,
    Parallelepiped,
    distance,
    generate_sites_in_region,
    # bounded lattice
    BoundedLattice,
    parallelepiped_region,
    dimension,
    # interact
    two_body_interaction_matrix,
    rydberg_interaction_matrix
    


    
include("lattice.jl")
include("region.jl")
include("bounded_lattice.jl")
include("interact.jl")
include("neighbors.jl")
include("visualize.jl")

end
