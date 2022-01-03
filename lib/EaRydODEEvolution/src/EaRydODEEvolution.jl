# Copyright 2020 QuEra Computing Inc. All rights reserved.

module EaRydODEEvolution

using Reexport
using Yao
using Adapt
using SparseArrays
using LinearAlgebra
using Configurations
using DiffEqCallbacks
using EaRydKrylovEvolution: AbstractTerm, AbstractSpace, EmulationOptions, storage_size, MemoryLayout, RealLayout, ComplexLayout
using OrdinaryDiffEq: OrdinaryDiffEq, Vern8, ODEProblem

@reexport using EaRydKrylovEvolution
export ShordingerEquation, ODEEvolution

struct EquationCache{H, Layout, S}
    hamiltonian::H
    layout::Layout
    state::S
end

function EquationCache(H::SparseMatrixCSC{Tv}, layout::ComplexLayout) where {Tv}
    state = Vector{Complex{real(Tv)}}(undef, size(H, 1))
    return EquationCache(H, layout, state)
end

function EquationCache(H::SparseMatrixCSC{Tv}, layout::RealLayout) where {Tv}
    state = Matrix{real(Tv)}(undef, size(H, 1), 2)
    return EquationCache(H, layout, state)
end

EquationCache(H::SparseMatrixCSC) = EquationCache(H, ComplexLayout())

struct ShordingerEquation{L, HTerm, Space, Cache <: EquationCache{<:Any, L}}
    layout::L
    hamiltonian::HTerm
    space::Space
    cache::Cache
end

function ShordingerEquation(h::AbstractTerm, space::AbstractSpace, cache::EquationCache)
    ShordingerEquation(cache.layout, h, space, cache)
end

Adapt.@adapt_structure ShordingerEquation
Adapt.@adapt_structure EquationCache

EaRydKrylovEvolution.storage_size(S::EquationCache) = storage_size(S.hamiltonian) + storage_size(S.state)

function Base.show(io::IO, m::MIME"text/plain", eq::ShordingerEquation)
    indent = get(io, :indent, 0)
    println(io, " "^indent, "Shordinger Equation:")
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

function (eq::ShordingerEquation)(dstate, state, p, t::Number) where L
    update_term!(eq.cache.hamiltonian, eq.hamiltonian(t), eq.space)
    mul!(eq.cache.state, eq.cache.hamiltonian, state)
    # @. dstate = -im * eq.cache.state
    update_dstate!(dstate, eq.cache.state, eq.layout)
    return
end

function update_dstate!(dstate::AbstractVector, state::AbstractVector, ::ComplexLayout)
    broadcast!(x->-im*x, dstate, state)
    return dstate
end

# real storage
# -im * (x + im*y)
# -im * x + y
# (y - x * im)
function update_dstate!(dstate::Matrix{<:Real}, state::Matrix{<:Real}, ::RealLayout)
    # real
    @inbounds for i in axes(state, 1)
        dstate[i, 1] = state[i, 2]
    end

    # imag
    @inbounds for i in axes(state, 1)
        dstate[i, 2] = -state[i, 1]
    end
    return dstate
end

function norm_preserve(resid, state, p, t)
    fill!(resid, 0)
    resid[1] = norm(state) - 1
    return
end

struct PieceWiseLinear{T}
    xs::Vector{T}
    ys::Vector{T}
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
struct ODEEvolution{P, Reg <: AbstractRegister, Eq <: ShordingerEquation, Prob <: ODEProblem, Options <: ODEOptions}
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

ODEEvolution(r::AbstractRegister, t, h::AbstractTerm; kw...) =
    ODEEvolution{real(Yao.datatype(r))}(r, t, h; kw...)

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
    layout = EaRydKrylovEvolution.MemoryLayout(r)
    if layout isa RealLayout
        isreal(h) || error("cannot use RealLayout for non-real hamiltonian")
    end
    options = ODEOptions(;kw...)
    # we do not convert the tspan to P since
    # the parameter function can relay on higher precision
    # and the performance of that usually doesn't matter
    start = EaRydKrylovEvolution.default_unit(μs, start)
    stop = EaRydKrylovEvolution.default_unit(μs, stop)
    time = (start, stop)
    reg = adapt(EaRydKrylovEvolution.PrecisionAdaptor(P), r)
    space = EaRydKrylovEvolution.get_space(r)

    # allocate cache
    # NOTE: on CPU we can do mixed type spmv
    # thus we use the smallest type we can get
    T = isreal(h) ? P : Complex{P}
    H = SparseMatrixCSC{T, Cint}(h(start+sqrt(eps(P))), space)
    cache = EquationCache(H, layout)
    eq = ShordingerEquation(h, space, cache)

    ode_prob = ODEProblem(
        eq, Yao.statevec(reg), time;
        save_everystep=false, save_start=false, alias_u0=true,
        progress=options.progress,
        progress_name=options.progress_name,
        progress_steps=options.progress_steps,
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
    printstyled(io, Base.format_bytes(EaRydKrylovEvolution.storage_size(prob.reg)); color=:yellow)
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

function EaRydKrylovEvolution.emulate!(prob::ODEEvolution)
    for _ in prob; end
    return prob
end

end
