module BloqadeExpr

using SparseArrays
using LinearAlgebra
using Adapt
using YaoAPI
using YaoBlocks
using LuxurySparse
using SparseMatricesCSR
using Base.Threads: nthreads
using MLStyle
using BitBasis
using LaTeXStrings
using ForwardDiff
using Unitful: Quantity, NoUnits, MHz, µm, uconvert
using InteractiveUtils: subtypes
using Base.Cartesian: @nexprs
using YaoBlocks: ChainBlock, PutBlock, TrivialGate, Subroutine, Scale, Daggered, Add, ControlBlock, TimeEvolution
using BloqadeLattices: BoundedLattice, rydberg_interaction_matrix


include("Lowlevel/Lowlevel.jl")

using .Lowlevel: Hamiltonian, SumOfLinop, ThreadedMatrix, storage_size, to_matrix, precision_type, highest_type, add_I, RegularLinop, isskewhermitian


export rydberg_h,
    rydberg_h_3,
    FullSpace,
    Subspace,
    fullspace,
    RydInteract,
    SumOfX, SumOfX_01, SumOfX_1r,
    SumOfXPhase, SumOfXPhase_01, SumOfXPhase_1r,
    SumOfZ, SumOfZ_01, SumOfZ_1r,
    SumOfN, SumOfN_1, SumOfN_r,
    XPhase, XPhase_01, XPhase_1r,
    PdPhase, PdPhase_01, PdPhase_1r,
    PuPhase, PuPhase_01, PuPhase_1r,
    X_01, X_1r,
    N_1, N_r,
    Pu_01, Pu_1r,
    Pd_01, Pd_1r,
    RydbergHamiltonian, RydbergHamiltonian3,
    get_rydberg_params,
    Op,
    attime,
    matrix_to_positions,
    storage_size,
    emulate!,
    precision_type,
    highest_type,

    to_matrix,
    add_I,
    RegularLinop,  # abstype
    SkewHermitian, # abstype
    isskewhermitian


include("assert.jl")
include("space.jl")
include("types.jl")
include("printings.jl")
include("mat.jl")


include("lower.jl")
include("units.jl")


include("interface.jl")
include("atoms.jl")


end
