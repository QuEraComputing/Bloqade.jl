


struct BloqadeDPSolver{Reg <: AbstractRegister, T, StateType, F, DPSolverType <: AbstractDPSolver{T, StateType, F}} <: AbstractDPSolver{T, StateType, F}
    tend::T
    reg::Reg
    dp_solver::DPSolverType
    function BloqadeDPSolver(
        tend::T,
        reg,
        dp_solver::AbstractDPSolver{T, StateType, F}
    ) where {T, StateType, F}
        new{typeof(reg), T, StateType, F, typeof(dp_solver)}(tend, reg, dp_solver)
    end
end

promote_tspan(tspan) = (zero(tspan), tspan)
promote_tspan(tspan::Tuple{A, B}) where {A, B} = promote_type(A, B).(tspan) 


function BloqadeDPSolver(reg::AbstractRegister, tspan, expr; solver_type=DP8Solver, copy_init=true,  kw...)
    nqudits(reg) == nqudits(expr) || throw(ArgumentError("number of qubits/sites does not match!"))
    reg = copy_init ? copy(reg) : reg
    tspan = promote_tspan(tspan)

    state = statevec(reg)
    space = YaoSubspaceArrayReg.space(reg)

    T = real(eltype(state))
    T = isreal(expr) ? T : Complex{T}
    eq = SchrodingerEquation(expr, Hamiltonian(T, expr, space))

    tspan_type = promote_type(real(eltype(state)), eltype(tspan))
    tstart, tend = tspan_type.(tspan) # promote tspan to T so Dual number works

    solver = solver_type(eq, tstart, state; kw...)
    BloqadeDPSolver(tend, reg, solver)
end


struct SchrodingerEquation{ExprType,H<:Hamiltonian}
    expr::ExprType
    hamiltonian::H
end


function (eq::SchrodingerEquation)(t::Real, state, dstate)
    fill!(dstate, zero(eltype(dstate)))
    for (f, term) in zip(eq.hamiltonian.fs, eq.hamiltonian.ts)
        mul!(dstate, term, state, -im * f(t), one(t))
    end
    return
end
