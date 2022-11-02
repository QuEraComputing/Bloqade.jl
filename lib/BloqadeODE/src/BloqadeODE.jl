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
using OrdinaryDiffEq: AutoSwitchCache, CompositeInterpolationData
using LinearAlgebra

using OrdinaryDiffEq:
    @logmsg,
    isadaptive,
    gamma_default,
    qmin_default,
    qmax_default,
    qsteady_min_default,
    qsteady_max_default,
    ODE_DEFAULT_NORM,
    ODE_DEFAULT_ISOUTOFDOMAIN,
    ODE_DEFAULT_UNSTABLE_CHECK,
    ODE_DEFAULT_PROG_MESSAGE,
    alg_extrapolates,
    DefaultInit,
    is_mass_matrix_alg,
    LogLevel,
    OrdinaryDiffEqAdaptiveAlgorithm,
    DAEAlgorithm,
    recursive_bottom_eltype,
    recursive_unitless_bottom_eltype,
    recursive_unitless_eltype,
    recursivecopy,
    initialize_tstops,
    initialize_saveat,
    initialize_d_discontinuities,
    uses_uprev,
    default_controller,
    DEOptions,
    OrdinaryDiffEqCompositeAlgorithm,
    InterpolationData,
    isdtchangeable,
    fsal_typeof,
    ODEIntegrator,
    initialize_callbacks!,
    handle_dt!

export SchrodingerProblem, SchrodingerEquation

include("problem.jl")
include("integrator_init.jl")

end
