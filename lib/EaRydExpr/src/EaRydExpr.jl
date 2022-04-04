module EaRydExpr

using SparseArrays
using LinearAlgebra
using YaoAPI
using YaoBlocks
using LuxurySparse
using MLStyle
using BitBasis
using LaTeXStrings
using Unitful: Quantity, NoUnits, MHz, µm, uconvert
using InteractiveUtils: subtypes
using Base.Cartesian: @nexprs
using YaoBlocks: ChainBlock, PutBlock, TrivialGate, Subroutine, Scale, Daggered, Add, ControlBlock, TimeEvolution

export rydberg_h, FullSpace, Subspace, fullspace,
    RydInteract, SumOfX, SumOfXPhase, SumOfZ, SumOfN, XPhase, Op,
    attime, matrix_to_positions, storage_size, emulate!

include("assert.jl")
include("space.jl")
include("types.jl")
include("printings.jl")
include("mat.jl")
include("linalg.jl")
include("lower.jl")
include("units.jl")
include("interface.jl")
include("atoms.jl")

end
