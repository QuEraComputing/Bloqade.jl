module BloqadeKrylov

using Adapt
using BloqadeExpr
using LinearAlgebra
using Configurations
using YaoArrayRegister
using YaoSubspaceArrayReg
using ForwardDiff
using BloqadeExpr: Hamiltonian, SumOfLinop
using ExponentialUtilities
using ProgressLogging
using GaussQuadrature

export emulate!, emulate_step!
export KrylovEvolution, Magnus4Evolution

export CFETEvolution, ACFETEvolution
export CFET2_1, CFET4_2, CFET6_5, CFET8_11

# utils.jl introduce a new type of Hamiltonian called ValHamiltonian:
# BloqadeExpr.Hamiltonian 
#   |-> BloqadeExpr.StepHamiltonian
#         |-> BloqadeKrylov.ValHamiltonian (new)
#include("utils.jl")

## following are expm_multiply
export onenormest, expm_multiply!, expmv!
include("onenormest.jl")
include("expm_multiply.jl")
include("expmv.jl")

include("common.jl")
include("krylov.jl")
include("magnus.jl")

## following are CFET
include("cfet.jl")
include("tables/cfet_tbl.jl")


## adaptive CFET 
include("adapt_cfet.jl")


if VERSION < v"1.7"
    include("patch.jl")
end



end
