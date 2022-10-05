# Copyright 2021 QuEra Computing Inc. All rights reserved.

module BloqadeLattices

using NearestNeighbors
using StatsBase
using LuxorGraphPlot
using LuxorGraphPlot: Point
using LuxorGraphPlot.Luxor: Colors

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
    distance

include("lattice.jl")
include("neighbors.jl")
include("visualize.jl")
include("region.jl")

end
