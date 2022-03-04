struct SchrodingerEquation{H, TC}
    hamiltonian::H
    term_cache::TC
end

function (eq::SchrodingerEquation)(dstate, state, p, t::Number) where L
    fill!(dstate, zero(eltype(dstate)))
    fs, hs = eq.term_cache.fs, eq.term_cache.hs
    for (f, h) in zip(fs, hs)
        # NOTE: currently we can expect all h
        # are preallocated constant matrices
        mul!(dstate, h, state, -im * f(t), one(t))
    end
    # NOTE: RealLayout is not supported
    # we will make it work automatically
    # later by using StructArrays
    return
end

function Base.show(io::IO, mime::MIME"text/plain", eq::SchrodingerEquation)
    indent = get(io, :indent, 0)
    tab(indent) = " "^indent
    print(io, tab(indent), "storage size: ")
    printstyled(io, Base.format_bytes(EaRydCore.storage_size(eq.term_cache)); color=:yellow)
    println(io)
    println(io, tab(indent), "expression:")
    EaRydCore.print_term(IOContext(io, :indent=>indent + 2), eq.hamiltonian)
    println(io)
    println(io)
end

struct SchrodingerProblem{Reg, EquationType <: ODEFunction, uType, tType, Kwargs} <: SciMLBase.AbstractODEProblem{uType, tType, true}
    reg::Reg
    f::EquationType
    state::uType # alias of reg.state
    u0::uType # a copy of input state
    tspan::tType

    # constant solver parameters
    kwargs::Kwargs

    p::Nothing # well make DiffEq happy

    function SchrodingerProblem(
        reg::AbstractRegister, tspan,
        hamiltonian::EaRydCore.AbstractTerm; kw...)

        nqubits(reg) == EaRydCore.nsites(hamiltonian) || throw(ArgumentError("number of qubits/sites does not match!"))
        # remove this after ArrayReg start using AbstractVector
        state = statevec(reg)
        space = EaRydCore.get_space(reg)
        tspan = SciMLBase.promote_tspan(tspan)
        # create term cache
        tc = EaRydCore.split_const_term(eltype(state), hamiltonian, space)

        eq = SchrodingerEquation(hamiltonian, tc)
        ode_f = ODEFunction(eq)

        default_ode_options = (
            save_everystep=false, save_start=false,
            save_on=false, dense=false,
        )
        kw = pairs(merge(default_ode_options, kw))

        new{typeof(reg), typeof(ode_f), typeof(state), typeof(tspan), typeof(kw)}(
            reg, ode_f, state, copy(state), tspan, kw, nothing
        )
    end
end

# multi-line printing
function Base.show(io::IO, mime::MIME"text/plain", prob::SchrodingerProblem)
    indent = get(io, :indent, 0)
    tab(indent) = " "^indent

    println(io, tab(indent), "SchrodingerProblem:")
    # state info
    println(io, tab(indent+2), "register info:")
    print(io, tab(indent+4), "type: ")
    printstyled(io, typeof(prob.reg); color=:green)
    println(io)

    print(io, tab(indent+4), "storage size: ")
    printstyled(io, Base.format_bytes(EaRydCore.storage_size(prob.reg)); color=:yellow)
    println(io)
    println(io)


    # tspan info
    println(io, tab(indent+2), "time span (μs): ", prob.tspan)
    println(io)

    # equation info
    println(io, tab(indent+2), "equation: ")
    show(IOContext(io, :indent=>indent+4), mime, prob.f.f)

    # options
    isempty(prob.kwargs) || println(io, tab(indent+2), "options:")
    order = [:algo, :dt, :adaptive, :progress]
    for key in order
        if haskey(prob.kwargs, key)
            println(io, tab(indent+4), key, ": ", repr(prob.kwargs[key]))
        end
    end

    for (key, value) in prob.kwargs
        key in order || println(io, tab(indent+4), key, ": ", repr(value))
    end
    # println(io)
end

function DiffEqBase.solve(prob::SchrodingerProblem, args...; sensealg=nothing, initial_state=nothing, kw...)
    if sensealg === nothing && haskey(prob.kwargs,:sensealg)
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

function EaRydCore.emulate!(prob::SchrodingerProblem)
    solve(prob, get(prob.kwargs, :algo, Vern8()))
    return prob
end
