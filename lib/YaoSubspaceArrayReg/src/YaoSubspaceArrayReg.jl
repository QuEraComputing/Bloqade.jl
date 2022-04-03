module YaoSubspaceArrayReg

using YaoAPI
using Random
using YaoArrayRegister
using EaRydExpr
using EaRydExpr: AbstractSpace
using BitBasis
using LinearAlgebra
using StatsBase
using YaoBlocks
using Base.Cartesian: @nexprs

export Subspace, SubspaceArrayReg, set_zero_state!,
    zero_state, rand_state, product_state,
    @bit_str, state, statevec, relaxedvec, isnormalized,
    nactive, nqubits, measure, measure!, space,
    fullspace

include("type.jl")
include("measure.jl")
include("instruct.jl")

end
