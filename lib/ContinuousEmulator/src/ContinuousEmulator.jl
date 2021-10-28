# Copyright 2020 QuEra Computing Inc. All rights reserved.

module ContinuousEmulator

using Reexport
using Yao
using Adapt
using SparseArrays
using LinearAlgebra
using Configurations
using DiffEqCallbacks
using RydbergEmulator: AbstractTerm, AbstractSpace, EmulationOptions, storage_size, MemoryLayout, RealLayout, ComplexLayout
using OrdinaryDiffEq: OrdinaryDiffEq, Vern8, ODEProblem

@reexport using RydbergEmulator
export ShordingerEquation, ContinuousEvolution

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

RydbergEmulator.storage_size(S::EquationCache) = storage_size(S.hamiltonian) + storage_size(S.state)

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

@option struct ContinuousOptions{Algo <: OrdinaryDiffEq.OrdinaryDiffEqAlgorithm} <: EmulationOptions
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
    ContinuousEvolution{P <: AbstractFloat}

Problem type for hamiltonian with time dependent parameters.
"""
struct ContinuousEvolution{P <: AbstractFloat, Reg <: AbstractRegister, Eq <: ShordingerEquation, Prob <: ODEProblem, Options <: ContinuousOptions}
    reg::Reg
    time::NTuple{2, P}
    eq::Eq
    ode_prob::Prob
    options::Options

    function ContinuousEvolution{P}(reg, time, eq, ode_prob, options) where P
        new{P, typeof(reg), typeof(eq), typeof(ode_prob), typeof(options)}(
            reg, time, eq, ode_prob, options
        )
    end
end

ContinuousEvolution(r::AbstractRegister, t, h::AbstractTerm; kw...) =
    ContinuousEvolution{real(Yao.datatype(r))}(r, t, h; kw...)

"""
    ContinuousEvolution{P}(r::AbstractRegister, t::Real, h::AbstractTerm; kw...)

Run the evolution for `t` μs, start from clock 0 μs, shorthand for

```julia
ContinuousEvolution{P}(r, (0, t), h; kw...)
```
"""
function ContinuousEvolution{P}(r::AbstractRegister, t::Real, h::AbstractTerm; kw...) where {P <: AbstractFloat}
    return ContinuousEvolution{P}(r, (zero(t), t), h; kw...)
end

"""
    ContinuousEvolution{P}(r::AbstractRegister, (start, stop), h::AbstractTerm; kw...) where {P <: AbstractFloat}

Create a `ContinuousEvolution` that defines the evolution of a hamiltonian `h` with time dependent parameters
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
function ContinuousEvolution{P}(r::AbstractRegister, (start, stop)::Tuple{<:Real, <:Real}, h::AbstractTerm; kw...) where {P <: AbstractFloat}
    layout = RydbergEmulator.MemoryLayout(r)
    if layout isa RealLayout
        isreal(h) || error("cannot use RealLayout for non-real hamiltonian")
    end
    options = ContinuousOptions(;kw...)
    # we do not convert the tspan to P since
    # the parameter function can relay on higher precision
    # and the performance of that usually doesn't matter
    start = RydbergEmulator.default_unit(μs, start)
    stop = RydbergEmulator.default_unit(μs, stop)
    time = (start, stop)
    reg = adapt(RydbergEmulator.PrecisionAdaptor(P), r)
    space = RydbergEmulator.get_space(r)

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
    return ContinuousEvolution{P}(reg, time, eq, ode_prob, options)
end

function Base.show(io::IO, mime::MIME"text/plain", prob::ContinuousEvolution{P}) where P
    indent = get(io, :indent, 0)
    tab(indent) = " "^indent
    println(io, tab(indent), "ContinuousEvolution{", P, "}:")
    # state info
    print(io, tab(indent), "  reg: ")
    printstyled(io, typeof(prob.reg); color=:green)
    println(io)
    print(io, tab(indent), "  reg storage: ")
    printstyled(io, Base.format_bytes(RydbergEmulator.storage_size(prob.reg)); color=:yellow)
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
        name = fieldname(ContinuousOptions, idx)
        print(io, tab(indent), "    $name: ", repr(getfield(prob.options, idx)))
        if idx != nfs
            println(io)
        end
    end
end

function RydbergEmulator.emulate!(prob::ContinuousEvolution)
    integrator = OrdinaryDiffEq.init(prob.ode_prob, prob.options.algo; reltol=prob.options.reltol, abstol=prob.options.abstol)
    for (idx, (u, t)) in enumerate(OrdinaryDiffEq.tuples(integrator))
        if iszero(mod(idx, prob.options.normalize_steps))
            normalize!(u)
        end
    end

    if prob.options.normalize_finally
        normalize!(prob.reg)
    end
    return prob
end

"""
    emulate!(r, t::Real, h::AbstractTerm; kw...)

Emulate the hamiltonian `h` for time duration `t`
on register `r` using ODE solvers. This is usually
the best way to emulate continuous hamiltonian
parameters.

# Arguments

- `r`: the register that contains the information about our quantum state of the system.
- `t`: the time duration to emulate.
- `h`: the hamiltonian terms.

# Selected Keyword Arguments

The keyword arguments are the same as DiffEq's ODEProblem interface,
we here only describe several commonly used keywords.

- `algo`: algorithm to use, see [Algorithm](#Algorithm) section for more discussion.
- `force_normalize`: whether to force the return statevector to be normalized, this
    is because for long-time emulation non-geometric methods will result in norm not
    equal to `1`, sometimes this error is not large enough to use a geometric method
    (which is usually slower), thus with some error we can normalize the state in the end
    using this option.
- `reltol`: relative tolerance, default is `1e-8`.
- `abstol`: absolute tolerance, default is `1e-8`.
- `progress`: whether showing progress bar, default is `false`.
- `progress_steps`: how much iteration steps should the progress estimate on,
    must be specified when `progress=true`.

# Algorithm

The default algorithm we usually use is the `Vern8` method, which provides
relative precise emulation for our Rydberg system. However, one may find lower
order approximation such as `Vern6` or `Vern7` can also be good for short time
emulation or smaller system. Or more general RK method such as `Tsit5` is also
a good candidate.

# Norm-Preservation via `ManifoldProjection`

One can enable norm-preservation via `ManifoldProjection`, the corresponding
callback function is implemented as `norm_preserve`.
"""
function RydbergEmulator.emulate!(r::Yao.AbstractRegister, t::Real, h::AbstractTerm; kw...)
    prob = ContinuousEvolution(r, t, h; kw...)
    emulate!(prob)
    return prob.reg
end

end
