"""
    struct SchrodingerEquation

Type for Schrodinger equation. A `SchrodingerEquation`
object is a callable object that has method `f(dstate, state, p, t)`
that fits into a standard ODE problem.
"""
struct SchrodingerEquation{ExprType,H<:Hamiltonian}
    expr::ExprType
    hamiltonian::H
end

Adapt.@adapt_structure SchrodingerEquation

function (eq::SchrodingerEquation)(dstate, state, p, t::Number) where {L}
    fill!(dstate, zero(eltype(dstate)))
    for (f, term) in zip(eq.hamiltonian.fs, eq.hamiltonian.ts)
        mul!(dstate, term, state, -im * f(t), one(t))
    end
    return
end

function Base.show(io::IO, mime::MIME"text/plain", eq::SchrodingerEquation)
    indent = get(io, :indent, 0)
    tab(indent) = " "^indent
    print(io, tab(indent), "storage size: ")
    printstyled(io, Base.format_bytes(storage_size(eq.hamiltonian)); color = :yellow)
    println(io)
    println(io, tab(indent), "expression:")
    show(IOContext(io, :indent => indent + 2), eq.expr)
    println(io)
    return println(io)
end

"""
    struct SchrodingerProblem
    SchrodingerProblem(reg, tspan, hamiltonian; kw...)

Define a Schrodinger equation problem that uses ODE solver from `OrdinaryDiffEq`
to solve the dynamics.

# Arguments

- `register`: required, the evolution problem register, can be a [`SubspaceArrayReg`](@ref) or an `ArrayReg`
    from `Yao`.
- `tspan`: required, a `(start, stop)` tuple or a single number `t`, the single value form `t` is equivalent
    to `(zero(t), t)`.
- `hamiltonian`: required, the evolution hamiltonian, can be created via [`rydberg_h`](@ref).

# Common Keyword Arguments

- `algo`: optional, algorithm to use, this only works for the `emulate!` interface.
    for `solve` or integrator interface, one will need to specify the algorithm explicitly.
- `progress`: print progress bar or not, this may effect the performance when problem scale is small, default is `true`.
- `progress_steps`: steps to update the progress bar, default is `5`.
- `reltol`: relative tolerance, default is 1e-8.
- `abstol`: absolute tolerance, default is 1e-8.

# Further References

For more ODE options, please refer to [Common Solver Options](https://diffeq.sciml.ai/stable/basics/common_solver_opts/).
The `SchrodingerProblem` type supports most of the standard DiffEq problem interface.
"""
struct SchrodingerProblem{Reg,EquationType<:ODEFunction,uType,tType,Algo,Kwargs} <:
       SciMLBase.AbstractODEProblem{uType,tType,true}
    reg::Reg
    f::EquationType
    state::uType # alias of reg.state
    u0::uType # a copy of input state
    tspan::tType
    algo::Algo

    # constant solver parameters
    kwargs::Kwargs

    p::Nothing # well make DiffEq happy
end

function SchrodingerProblem(reg::AbstractRegister, tspan, expr; algo=AutoVern9(Rodas5(autodiff=false)), kw...)
    nqudits(reg) == nqudits(expr) || throw(ArgumentError("number of qubits/sites does not match!"))
    # remove this after ArrayReg start using AbstractVector
    state = statevec(reg)
    space = YaoSubspaceArrayReg.space(reg)
    tspan = SciMLBase.promote_tspan(tspan)
    # create term cache
    # always follow register element-type
    T = real(eltype(state))
    T = isreal(expr) ? T : Complex{T}
    eq = SchrodingerEquation(expr, Hamiltonian(T, expr, space))
    ode_f = ODEFunction(eq)

    tspan_type = promote_type(real(eltype(state)), eltype(tspan))
    tspan = tspan_type.(tspan) # promote tspan to T so Dual number works

    default_ode_options = (
        save_everystep = false,
        save_start = false,
        save_on = false,
        dense = false,
        reltol=1e-8,
        abstol=1e-8,
    )
    kw = pairs(merge(default_ode_options, kw))

    return SchrodingerProblem{typeof(reg),typeof(ode_f),typeof(state),typeof(tspan),typeof(algo), typeof(kw)}(
        reg,
        ode_f,
        state,
        copy(state),
        tspan,
        algo,
        kw,
        nothing,
    )
end

function Adapt.adapt_structure(to, x::SchrodingerProblem)
    d_reg = adapt(to, x.reg)
    state = statevec(d_reg)
    return SchrodingerProblem(
        d_reg, adapt(to, x.f), state,
        copy(state), x.tspan, x.algo, x.kwargs, x.p
    )
end

# multi-line printing
function Base.show(io::IO, mime::MIME"text/plain", prob::SchrodingerProblem)
    indent = get(io, :indent, 0)
    tab(indent) = " "^indent

    println(io, tab(indent), "SchrodingerProblem:")
    # state info
    println(io, tab(indent + 2), "register info:")
    print(io, tab(indent + 4), "type: ")
    printstyled(io, typeof(prob.reg); color = :green)
    println(io)

    print(io, tab(indent + 4), "storage size: ")
    printstyled(io, Base.format_bytes(storage_size(prob.reg)); color = :yellow)
    println(io)
    println(io)

    # tspan info
    println(io, tab(indent + 2), "time span (μs): ", prob.tspan)
    println(io)

    # equation info
    println(io, tab(indent + 2), "equation: ")
    show(IOContext(io, :indent => indent + 4), mime, prob.f.f)
    println(io, tab(indent + 4), "algorithm: ", repr(prob.algo))

    # options
    isempty(prob.kwargs) || println(io, tab(indent + 2), "options:")
    order = [:dt, :adaptive, :progress]
    for key in order
        if haskey(prob.kwargs, key)
            println(io, tab(indent + 4), key, ": ", repr(prob.kwargs[key]))
        end
    end

    for (key, value) in prob.kwargs
        key in order || println(io, tab(indent + 4), key, ": ", repr(value))
    end
    # println(io)
end

function DiffEqBase.solve(prob::SchrodingerProblem, args...; sensealg = nothing, initial_state = nothing, kw...)
    if sensealg === nothing && haskey(prob.kwargs, :sensealg)
        sensealg = prob.kwargs[:sensealg]
    end

    # update initial state
    if initial_state !== nothing
        initial_state isa AbstractRegister ||
            throw(ArgumentError("initial_state must be a register, got $(typeof(initial_state))"))
        u0 = statevec(initial_state)
    else
        u0 = prob.u0
    end
    return DiffEqBase.solve_up(prob, sensealg, u0, nothing, args...; kw...)
end

DiffEqBase.get_concrete_problem(prob::SchrodingerProblem, isadapt; kw...) = prob

function BloqadeExpr.emulate!(prob::SchrodingerProblem)
    solve(prob, prob.algo)
    return prob
end
