module BloqadeDormandPrince

using BloqadeExpr: BloqadeExpr, Hamiltonian, nqudits
using YaoArrayRegister: AbstractRegister, ArrayReg, statevec
using YaoSubspaceArrayReg: YaoSubspaceArrayReg
using Reexport: @reexport
@reexport using DormandPrince
using DormandPrince: integrate_core!
using LinearAlgebra: mul!

export SchrodingerProblem

include("types.jl")
include("impl.jl")

end # BloqadeDormandPrince
