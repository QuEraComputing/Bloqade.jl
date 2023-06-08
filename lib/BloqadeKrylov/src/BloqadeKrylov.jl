module BloqadeKrylov

using Adapt
using BloqadeExpr
using LinearAlgebra
using Configurations
using YaoArrayRegister
using YaoSubspaceArrayReg
using BloqadeExpr: Hamiltonian, StepHamiltonian
using ExponentialUtilities
using ProgressLogging

export emulate!, emulate_step!
export KrylovEvolution, Magnus4Evolution

export CFETEvolution
export CFET42

# utils.jl introduce a new type of Hamiltonian called ValHamiltonian:
# BloqadeExpr.Hamiltonian 
#   |-> BloqadeExpr.StepHamiltonian
#         |-> BloqadeKrylov.ValHamiltonian
include("utils.jl")


include("expmv.jl")
include("common.jl")
include("krylov.jl")
include("magnus.jl")

## following are CFET
include("cfet.jl")
include("tables/cfet42_tbl.jl")


if VERSION < v"1.7"
    include("patch.jl")
end



end
