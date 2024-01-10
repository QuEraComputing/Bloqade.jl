module BloqadeODE

using Adapt
using Reexport
using SciMLBase
using DiffEqBase
using YaoArrayRegister
using YaoSubspaceArrayReg
@reexport using BloqadeExpr
@reexport using OrdinaryDiffEq
using BloqadeExpr: Hamiltonian
using LinearAlgebra

export SchrodingerProblem, SchrodingerEquation

include("problem.jl")

end
