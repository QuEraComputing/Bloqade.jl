module RydbergEmulator

using Printf
using BitBasis
using SparseArrays
using LuxurySparse
using LightGraphs
using LinearAlgebra
using OrderedCollections
using ExponentialUtilities
using StatsBase
using Random
using Printf

import Yao
using Yao: measure

using LinearAlgebra: BlasReal, BlasComplex

export RydInteract, RydAtom, XTerm, ZTerm, Hamiltonian, EmulatorCache, RydbergReg
export to_matrix!, update_term!, simple_rydberg, rydberg_h, rydatoms, rand_atoms,
    unit_disk_graph, rand_unit_disk_graph, emulate!, emulate, square_lattice, set_zero_state!, measure

include("utils.jl")
include("atoms.jl")
include("subspace.jl")
include("hamiltonian.jl")

include("register.jl")
include("measure.jl")

include("unit_disk.jl")
include("emulate.jl")
include("qaoa_mis.jl")

include("deprecations.jl")

end # module
