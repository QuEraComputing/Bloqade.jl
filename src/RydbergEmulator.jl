module RydbergEmulator

using Printf
using BitBasis
using SparseArrays
using LuxurySparse
using LightGraphs
using LinearAlgebra
using OrderedCollections
using ExponentialUtilities
using Configurations
using DelimitedFiles
using StatsBase
using Random
using Printf

import Yao
using Yao: measure, zero_state

using LinearAlgebra: BlasReal, BlasComplex

export RydInteract, RydAtom, XTerm, ZTerm, Hamiltonian, EmulatorCache, RydbergReg, Subspace
export to_matrix!, update_term!, simple_rydberg, rydberg_h, rydatoms, rand_atoms, read_atoms, write_atoms,
    unit_disk_graph, rand_unit_disk_graph, emulate!, emulate,
    square_lattice, set_zero_state!, blockade_subspace,
    # reexport from Yao
    measure, zero_state

include("utils.jl")
include("atoms.jl")
include("subspace.jl")
include("hamiltonian.jl")

include("register.jl")
include("measure.jl")

include("unit_disk.jl")
include("emulate.jl")
include("qaoa_mis.jl")

include("serialize.jl")

include("deprecations.jl")

end # module
