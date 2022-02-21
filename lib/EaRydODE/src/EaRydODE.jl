# Copyright 2020 QuEra Computing Inc. All rights reserved.

module EaRydODE

using Reexport
using Yao
using Adapt
using SparseArrays
using LinearAlgebra
using Configurations
using DiffEqCallbacks
using EaRydCore: AbstractTerm, AbstractSpace, EmulationOptions,
    storage_size, nsites, MemoryLayout, RealLayout, ComplexLayout,
    split_term
using OrdinaryDiffEq: OrdinaryDiffEq, ODEProblem

@reexport using EaRydCore
@reexport using OrdinaryDiffEq: Vern6, Vern7, Vern8, VCABM, AB3
export SchrodingerEquation, ODEEvolution

struct EquationCache{H, Layout, S}
    hamiltonian::H
    layout::Layout
    state::S
end

function EquationCache(::Type{Tv}, h::AbstractTerm, space::AbstractSpace, layout::ComplexLayout) where {Tv}
    state = Vector{Complex{real(Tv)}}(undef, size(H, 1))
    tc = split_term(Tv, h, space)
    return EquationCache(tc, layout, state)
end

function EquationCache(::Type{Tv}, h::AbstractTerm, space::AbstractSpace, layout::RealLayout) where {Tv}
    state = Matrix{real(Tv)}(undef, size(H, 1), 2)
    tc = split_term(Tv, h, space)
    return EquationCache(tc, layout, state)
end

struct SchrodingerEquation{L, HTerm, Space, Cache <: EquationCache{<:Any, L}}
    layout::L
    hamiltonian::HTerm
    space::Space
    cache::Cache
end

function SchrodingerEquation(h::AbstractTerm, space::AbstractSpace, cache::EquationCache)
    SchrodingerEquation(cache.layout, h, space, cache)
end

Adapt.@adapt_structure SchrodingerEquation
Adapt.@adapt_structure EquationCache

function EaRydCore.storage_size(S::EquationCache)
    return storage_size(S.hamiltonian) + storage_size(S.state)
end

function Base.show(io::IO, m::MIME"text/plain", eq::SchrodingerEquation)
    indent = get(io, :indent, 0)
    println(io, " "^indent, "Schrödinger Equation:")
    print(io, " "^indent, "  Storage Size: ")
    printstyled(io, Base.format_bytes(storage_size(eq.cache)); color=:green)
    println(io)
    print(io, " "^indent, "  State Storage: ")
    printstyled(io, typeof(eq.cache.state); color=:green)
    println(io)
    print(io, " "^indent, "  Hamiltonian Storage: ")
    printstyled(io, typeof(eq.cache.hamiltonian); color=:green)
    println(io)

    show(IOContext(io, :indent=>indent+2), m, eq.hamiltonian)
end

function (eq::SchrodingerEquation)(dstate, state, p, t::Number) where L
    fs, hs = eq.cache.hamiltonian.fs, eq.cache.hamiltonian.hs
    for (f, h) in zip(fs, hs)
        # NOTE: currently we can expect all h
        # are preallocated constant matrices
        mul!(dstate, h, state, f(t), one(t))
    end
    # NOTE: RealLayout is not supported
    # we will make it work automatically
    # later by using StructArrays
    lmul!(-im, dstate)
    return dstate
end

@option struct ODEOptions{Algo <: OrdinaryDiffEq.OrdinaryDiffEqAlgorithm} <: EmulationOptions
    algo::Algo = Vern8()
    progress::Bool = false
    progress_steps::Int = 5
    progress_name::String = "ODE"
    reltol::Float64 = 1e-8
    abstol::Float64 = 1e-8
    normalize_steps::Int = 5
    normalize_finally::Bool = true
end

"""
    ODEEvolution{P}

Problem type for hamiltonian with time dependent parameters.
"""
struct ODEEvolution{P, Reg <: AbstractRegister, Eq <: SchrodingerEquation, Prob <: ODEProblem, Options <: ODEOptions}
    reg::Reg
    time::NTuple{2, P}
    eq::Eq
    ode_prob::Prob
    options::Options

    function ODEEvolution{P}(reg, time, eq, ode_prob, options) where P
        new{P, typeof(reg), typeof(eq), typeof(ode_prob), typeof(options)}(
            reg, time, eq, ode_prob, options
        )
    end
end

function ODEEvolution(r::AbstractRegister, t, h::AbstractTerm; kw...)
    PR = real(Yao.datatype(r))
    PH = real(eltype(h(zero(t))))
    P = promote_type(PR, PH)
    return ODEEvolution{P}(r, t, h; kw...)
end

"""
    ODEEvolution{P}(r::AbstractRegister, t::Real, h::AbstractTerm; kw...)

Run the evolution for `t` μs, start from clock 0 μs, shorthand for

```julia
ODEEvolution{P}(r, (0, t), h; kw...)
```
"""
function ODEEvolution{P}(r::AbstractRegister, t::Real, h::AbstractTerm; kw...) where {P}
    return ODEEvolution{P}(r, (zero(t), t), h; kw...)
end

"""
    ODEEvolution{P}(r::AbstractRegister, (start, stop), h::AbstractTerm; kw...) where {P <: AbstractFloat}

Create a `ODEEvolution` that defines the evolution of a hamiltonian `h` with time dependent parameters
to evolve from `start` to `stop` using an ODE solver.

# Arguments

- `P`: optional, a type parameter that sets the problem precision type, default is
    the same as the `Yao.datatype` of given `register`.
- `register`: required, the evolution problem register, can be a [`RydbergReg`](@ref) or an `ArrayReg`
    from `Yao`.
- `(start, stop)`: required, the evolution interval.
- `h`: required, the evolution hamiltonian.

# Keyword Arguments

- `algo`: algorithm to use, default is `Vern8`, check DiffEq documentation for more details.
- `progress`: print progress bar or not, this may effect the performance when problem scale is small, default is `true`.
- `progress_steps`: steps to update the progress bar, default is `5`.
- `reltol`: relative tolerance, default is 1e-8.
- `abstol`: absolute tolerance, default is 1e-8.
- `normalize_steps`: steps to run normalization on the state, default is `5`.
"""
function ODEEvolution{P}(r::AbstractRegister, (start, stop)::Tuple{<:Real, <:Real}, h::AbstractTerm; kw...) where {P}
    nqubits(r) == nsites(h) || error("number of sites does not match")
    
    layout = EaRydCore.MemoryLayout(r)
    if layout isa RealLayout
        isreal(h) || error("cannot use RealLayout for non-real hamiltonian")
    end

    main_options, ode_options = [], []
    for (k, v) in kw
        if k in fieldnames(ODEOptions)
            push!(main_options, k=>v)
        else
            push!(ode_options, k=>v)
        end
    end

    options = ODEOptions(;main_options...)
    # we do not convert the tspan to P since
    # the parameter function can relay on higher precision
    # and the performance of that usually doesn't matter
    start = EaRydCore.default_unit(μs, start)
    stop = EaRydCore.default_unit(μs, stop)
    time = (start, stop)
    reg = adapt(EaRydCore.PrecisionAdaptor(P), r)
    space = EaRydCore.get_space(r)

    # allocate cache
    # NOTE: on CPU we can do mixed type spmv
    # thus we use the smallest type we can get
    T = isreal(h) ? P : Complex{P}
    H = SparseMatrixCSC{T, Cint}(h(start+sqrt(eps(P))), space)
    cache = EquationCache(H, layout)
    eq = SchrodingerEquation(h, space, cache)

    ode_prob = ODEProblem(
        eq, Yao.statevec(reg), time;
        save_everystep=false, save_start=false, alias_u0=true,
        options.progress,
        options.progress_name,
        options.progress_steps,
        ode_options...
    )
    return ODEEvolution{P}(reg, time, eq, ode_prob, options)
end

function Base.show(io::IO, mime::MIME"text/plain", prob::ODEEvolution{P}) where P
    indent = get(io, :indent, 0)
    tab(indent) = " "^indent
    println(io, tab(indent), "ODEEvolution{", P, "}:")
    # state info
    print(io, tab(indent), "  reg: ")
    printstyled(io, typeof(prob.reg); color=:green)
    println(io)
    print(io, tab(indent), "  reg storage: ")
    printstyled(io, Base.format_bytes(EaRydCore.storage_size(prob.reg)); color=:yellow)
    println(io)
    println(io)

    # time span
    println(io, tab(indent), "  timespan: ", prob.ode_prob.tspan)
    # equation
    println(io, tab(indent), "  equation: ")
    show(IOContext(io, :indent=>indent+4), mime, prob.eq)

    println(io)
    println(io)
    println(io, tab(indent), "  options: ")

    nfs = nfields(prob.options)
    for idx in 1:nfs
        name = fieldname(ODEOptions, idx)
        print(io, tab(indent), "    $name: ", repr(getfield(prob.options, idx)))
        if idx != nfs
            println(io)
        end
    end
end

function get_integrator(prob::ODEEvolution)
    return OrdinaryDiffEq.init(
        prob.ode_prob, prob.options.algo;
        reltol=prob.options.reltol, abstol=prob.options.abstol
    )
end

function emulate_step!(prob::ODEEvolution, ret)
    if prob.options.normalize_finally && ret === nothing
        normalize!(prob.reg)
    end

    ret === nothing && return
    step, reg =ret[1][1], prob.reg

    if iszero(mod(step, prob.options.normalize_steps))
        normalize!(reg)
    end
    return (;step, reg, clock=ret[1][2])
end

function Base.iterate(prob::ODEEvolution)
    integrator = enumerate(get_integrator(prob))
    ret = iterate(integrator)
    info = emulate_step!(prob, ret)
    return info, (integrator, ret[2])
end

function Base.iterate(prob::ODEEvolution, (integrator, st))
    ret = iterate(integrator, st)
    info = emulate_step!(prob, ret)
    ret === nothing && return
    return info, (integrator, ret[2])
end

function EaRydCore.emulate!(prob::ODEEvolution)
    for _ in prob; end
    return prob
end

end
