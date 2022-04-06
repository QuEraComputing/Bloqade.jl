module EaRydKrylov

using EaRydExpr
using LinearAlgebra
using Configurations
using YaoArrayRegister
using YaoSubspaceArrayReg
using EaRydExpr: Hamiltonian, StepHamiltonian
using ExponentialUtilities
using ProgressLogging

export KrylovEvolution, emulate!, emulate_step!

include("expmv.jl")
include("emulate.jl")

end
