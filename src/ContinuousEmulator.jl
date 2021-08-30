# Copyright 2020 QuEra Computing Inc. All rights reserved.

module ContinuousEmulator

using Reexport
@reexport using RydbergEmulator
using OrdinaryDiffEq
using Yao
using Adapt
using SparseArrays
using LinearAlgebra
using DiffEqCallbacks
using RydbergEmulator: AbstractTerm

export ShordingerEquation

struct EquationCache{H, S}
    hamiltonian::H
    state::S
end

# CPU
function EquationCache(H::SparseMatrixCSC{Tv}) where Tv
    state = Vector{Complex{real(Tv)}}(undef, size(H, 1))
    return EquationCache(H, state)
end

struct ShordingerEquation{S, HTerm, HMatrix, State}
    subspace::S
    hamiltonian::HTerm
    cache::EquationCache{HMatrix, State}
end

Adapt.@adapt_structure ShordingerEquation
Adapt.@adapt_structure EquationCache

# default is always CPU, we use cu(x) interface for CUDA via Adapt

"""
    ShordingerEquation([T=Float64, s::Subspace], h::AbstractTerm[, cache=cpu_equation_cache(T, h)])

Create a Shordinger equation for given hamiltonian expression `h`.

# Arguments

- `T <: Real`: eltype of the hamiltonian matrix, default is `Float64`
    when the hamiltonian is a real hermitian, or `ComplexF64` when the hamiltonian
    is a complex hermitian to acheive best performance.
- `s::Subspace`: subspace, can be used to specify an
    independent-set subspace for the simulation, optional.
- `h::AbstractTerm`: required, the hamiltonian.
- `cache::EquationCache`: the equation evaluation cache, default is `cpu_equation_cache(T, h)`.

# The Callable Method

A `ShordingerEquation` object has a callable method

    (eq::ShordingerEquation)(dstate, state, p, t)

Can be further used in an ODE solver such as `OrdinaryDiffEq`.
The form `(dstate, state, p, t)` is the inplace definition
of the equation, see also
[DiffEq:Solving-Systems-of-Equations](https://diffeq.sciml.ai/latest/tutorials/ode_example/#Example-2:-Solving-Systems-of-Equations)
"""
ShordingerEquation(h::AbstractTerm) = ShordingerEquation(decide_eltype(h), h)
ShordingerEquation(h::AbstractTerm, cache::EquationCache) = ShordingerEquation(decide_eltype(h), h, cache)
ShordingerEquation(s::Subspace, h::AbstractTerm) = ShordingerEquation(decide_eltype(h), s, h)

decide_eltype(h::AbstractTerm) = isreal(h) ? Float64 : ComplexF64

function ShordingerEquation(::Type{T}, h::AbstractTerm, cache::EquationCache=cpu_equation_cache(T, h)) where T
    ShordingerEquation(nothing, h, cache)
end

function ShordingerEquation(::Type{T}, s::Subspace, h::AbstractTerm, cache::EquationCache=cpu_equation_cache(T, h, s)) where T
    ShordingerEquation(s, h, cache)
end

# NOTE: we choose the matrix eltype automatically when it can be real
# since SparseMatrix{<:Real} * Vector{<:Complex} is always faster than
# SparseMatrix{<:Complex} * Vector{<:Complex}
function cpu_equation_cache(::Type{T}, h::AbstractTerm, s::Subspace) where T
    EquationCache(SparseMatrixCSC{T}(h(1e-2), s))
end

function cpu_equation_cache(::Type{T}, h::AbstractTerm) where T
    EquationCache(SparseMatrixCSC{T}(h(1e-2)))
end

estimate_size(S::SparseMatrixCSC) = sizeof(S.colptr) + sizeof(S.rowval) + sizeof(S.nzval)
estimate_size(S) = sizeof(S)
estimate_size(S::EquationCache) = estimate_size(S.hamiltonian) + estimate_size(S.state)

function Base.show(io::IO, m::MIME"text/plain", eq::ShordingerEquation)
    indent = get(io, :indent, 0)
    println(io, " "^indent, "Shordinger Equation:")
    print(io, " "^indent, "  Storage Size: ")
    printstyled(io, Base.format_bytes(estimate_size(eq.cache)); color=:green)
    println(io)
    print(io, " "^indent, "  State Storage: ")
    printstyled(io, typeof(eq.cache.state); color=:green)
    println(io)
    print(io, " "^indent, "  Hamiltonian Storage: ")
    printstyled(io, typeof(eq.cache.hamiltonian); color=:green)
    println(io)

    show(IOContext(io, :indent=>indent+2), m, eq.hamiltonian)
end

function (eq::ShordingerEquation{<:Subspace})(dstate, state, p, t::Number)
    update_term!(eq.cache.hamiltonian, eq.hamiltonian(t), eq.subspace)
    mul!(eq.cache.state, eq.cache.hamiltonian, state)
    @. dstate = -im * eq.cache.state
    return
end

function (eq::ShordingerEquation{Nothing})(dstate, state, p, t::Number)
    update_term!(eq.cache.hamiltonian, eq.hamiltonian(t))
    mul!(eq.cache.state, eq.cache.hamiltonian, state)
    @. dstate = -im * eq.cache.state
    return
end

function norm_preserve(resid, state, p, t)
    fill!(resid, 0)
    resid[1] = norm(state) - 1
    return
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
function RydbergEmulator.emulate!(r::Yao.ArrayReg, t::Real, h::AbstractTerm;
        algo=Vern8(), force_normalize=true, reltol=1e-8, abstol=1e-8,
        kwargs...
    )
    prob = ODEProblem(
        ShordingerEquation(h),
        vec(r.state), (zero(t), t);
        save_everystep=false, save_start=false, alias_u0=true, kwargs...
    )
    result = solve(prob, algo; reltol, abstol)

    # force normalize
    if force_normalize
        r.state ./= norm(vec(r.state))
    end
    return r
end

function RydbergEmulator.emulate!(r::RydbergReg, t::Real, h::AbstractTerm;
        algo=Vern8(), force_normalize=true, reltol=1e-8, abstol=1e-8,
        kwargs...
    )

    prob = ODEProblem(
        ShordingerEquation(r.subspace, h),
        vec(r.state), (zero(t), t);
        save_everystep=false, save_start=false, alias_u0=true, kwargs...
    )
    result = solve(prob, algo; reltol, abstol)

    # force normalize
    if force_normalize
        r.state ./= norm(vec(r.state))
    end
    return r
end

end
