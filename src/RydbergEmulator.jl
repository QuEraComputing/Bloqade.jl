module RydbergEmulator

using Random
using Printf
using Unitful
using BitBasis
using SparseArrays
using LuxurySparse
using LightGraphs
using LinearAlgebra
using OrderedCollections
using Configurations
using DelimitedFiles
using EliminateGraphs
using StatsBase

import Yao
using Yao: measure, zero_state
using Unitful: Quantity, uconvert, MHz, µm, μs, ns
using LinearAlgebra: BlasReal, BlasComplex

export RydInteract, RydAtom, XTerm, ZTerm, NTerm, Hamiltonian, EmulatorCache, RydbergReg, Subspace
export to_matrix!, update_term!, simple_rydberg, rydberg_h, rydatoms, rand_atoms, read_atoms, write_atoms,
    unit_disk_graph, rand_unit_disk_graph, emulate!, emulate,
    square_lattice, set_zero_state!, blockade_subspace, is_independent_set, to_independent_set!,
    # reexport from Yao
    measure, zero_state,
    # units
    MHz, μs, ns


# TODO: move to ExpV.jl
include("expv/exp.jl")
include("expv/expv.jl")
include("expv/krylov.jl")

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
include("mis.jl")
include("deprecations.jl")

end # module
