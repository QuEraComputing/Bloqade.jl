module EaRydCore

using Random
using Printf
using Unitful
using BitBasis
using UUIDs
using Adapt
using SparseArrays
using LuxurySparse
using Graphs
using LinearAlgebra
using OrderedCollections
using Configurations
using DelimitedFiles
using EliminateGraphs
using ProgressLogging
using Polyester
using StatsBase
using ThreadsX
using Transducers
using LaTeXStrings
using StaticArrays
using Reexport

using Transducers: OnInit
using MLStyle: @match, @switch
using LinearAlgebra: BlasInt, BlasReal, BlasComplex

@reexport using Unitful: Quantity, uconvert, MHz, µm, μs, ns
@reexport using Yao

export RydInteract, RydAtom, XTerm, ZTerm, NTerm,
    Hamiltonian, KrylovEmulationCache,
    RydbergReg, RealLayout, ComplexLayout, MemoryLayout,
    FullSpace, fullspace, Subspace, SubspaceMap,
    is_time_dependent,
    update_term!, simple_rydberg, rydberg_h, rydatoms,
    rand_atoms, read_atoms, write_atoms, read_subspace,
    write_subspace, unit_disk_graph, rand_unit_disk_graph,
    KrylovEvolution, emulate!, trotterize,
    mean_rydberg, count_vertices, mean, gibbs_loss, logsumexp,
    square_lattice, set_zero_state!, blockade_subspace, independent_set_subspace,
    is_independent_set, to_independent_set!, to_independent_set, add_vertices, add_vertices!,
    add_random_vertices, independent_set_probabilities,
    mis_postprocessing, Op,
    # observables
    rydberg_density

# NOTE: remove this after expv get fixed
include("expmv.jl")

# NOTE: remove this after BQCESubroutine beta version is released
include("bsubspace.jl")

include("utils.jl")
include("atoms.jl")
include("subspace.jl")
include("hamiltonian/hamiltonian.jl")

include("register.jl")
include("measure.jl")
include("instructs.jl")
include("mat.jl")

include("unit_disk.jl")
include("emulate/emulate.jl")

include("serialize.jl")
include("mis.jl")
include("deprecations.jl")

include("observables.jl")

include("print_latex.jl")


end # module
