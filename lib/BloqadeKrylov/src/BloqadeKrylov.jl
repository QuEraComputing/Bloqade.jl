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
export KrylovEvolution, Magnus4Evolution, CFET42Evolution


include("expmv.jl")
include("common.jl")
include("krylov.jl")
include("magnus.jl")
include("cfet.jl")


if VERSION < v"1.7"
    include("patch.jl")
end



end
