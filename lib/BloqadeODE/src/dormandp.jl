
struct BloqadeSolver{Reg <: AbstractRegister, T, StateType, F, DPSolverType <: AbstractDPSolver{T, StateType, F}} <: AbstractDPSolver{T, StateType, F}
    reg::Reg
    dp_solver::DPSolverType
    function BloqadeSolver(
        reg,
        dp_solver::AbstractDPSolver{T, StateType, F}
    ) where {T, StateType, F}
        new{typeof(reg), T, StateType, F, typeof(dp_solver)}(reg, dp_solver)
    end
end

function BloqadeSolver(reg::AbstractRegister, tstart::Real, expr; solver_type=DP8Solver, copy_init=true,  kw...)
    nqudits(reg) == nqudits(expr) || throw(ArgumentError("number of qubits/sites does not match!"))
    # remove this after ArrayReg start using AbstractVector
    reg = copy_init ? copy(reg) : reg

    state = statevec(reg)
    space = YaoSubspaceArrayReg.space(reg)

    T = real(eltype(state))
    T = isreal(expr) ? T : Complex{T}
    eq = SchrodingerEquation(expr, Hamiltonian(T, expr, space))

    solver = solver_type(eq, convert(T, tstart), state; kw...)
    BloqadeSolver(reg, solver)
end

DormandPrince.get_current_state(solver::BloqadeSolver) = solver.reg
DormandPrince.integrate_core!(solver::BloqadeSolver, time::Real) = integrate_core!(solver.dp_solver, time)
    
