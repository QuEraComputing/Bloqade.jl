module RydbergEmulator

using Random
using Printf
using Unitful
using BitBasis
using MLStyle
using UUIDs
using SparseArrays
using LuxurySparse
using LightGraphs
using LinearAlgebra
using OrderedCollections
using Configurations
using DelimitedFiles
using EliminateGraphs
using ProgressLogging
using StatsBase
using ThreadsX
using StaticArrays

import Yao
using Yao: measure, zero_state
using Unitful: Quantity, uconvert, MHz, µm, μs, ns
using LinearAlgebra: BlasReal, BlasComplex

export RydInteract, RydAtom, XTerm, ZTerm, NTerm, Hamiltonian, EmulatorCache, RydbergReg, Subspace,
    PulseJob, Pulse, HyperfinePulse, RydbergPulse
export to_matrix!, update_term!, simple_rydberg, rydberg_h, rydatoms, rand_atoms, read_atoms, write_atoms,
    unit_disk_graph, rand_unit_disk_graph, emulate!, emulate,
    square_lattice, set_zero_state!, blockade_subspace, is_independent_set, to_independent_set!,
    # reexport from Yao
    measure, zero_state, reduced_mean_rydberg,
    # units
    MHz, μs, ns

# NOTE: remove this after expv get fixed
include("expmv.jl")

# NOTE: remove this after BQCESubroutine beta version is released
include("bsubspace.jl")

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
include("schema.jl")
include("deprecations.jl")

end # module
