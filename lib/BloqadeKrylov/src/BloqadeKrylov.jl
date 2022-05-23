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

export KrylovEvolution, emulate!, emulate_step!

include("expmv.jl")
include("emulate.jl")

if VERSION < v"1.7"
    include("patch.jl")
end

end
