function emulate_routine!(r::RydbergReg, t::Number, h::AbstractTerm, cache::AbstractMatrix)
    st = vec(r.state)
    update_term!(cache, h, r.subspace)
    expmv!(-im * t, cache, st)
    return r
end

function emulate_routine!(r::Yao.ArrayReg, t::Number, h::AbstractTerm, cache::AbstractMatrix)
    st = vec(r.state)
    update_term!(cache, h)
    expmv!(-im * t, cache, st)
    return r
end

"""
    EmulatorCache

Cache type for the emulator.
"""
struct EmulatorCache{C}
    H::C
end

function EmulatorCache(::Type{T}, H::AbstractTerm) where {T <: Real}
    return EmulatorCache(SparseMatrixCSC{Complex{T}}(H))
end

function EmulatorCache(::Type{T}, H::AbstractTerm, s::Subspace) where {T <: Real}
    return EmulatorCache(SparseMatrixCSC{Complex{T}}(H, s))
end

EmulatorCache(H::AbstractTerm, xs...) = EmulatorCache(eltype(H), H, xs...)

"""
    emulate!(r, ts, hs[; cache=EmulatorCache(ts, hs)])

Emulate the time evolution of a sequence of Hamiltonians `hs` of time length `ts` for each Hamiltonian.
An optional argument can be feeded to preallocate the memory.

    emulate!(r, t::Real, h::AbstractTerm[; algo=Tsit5(), cache=<default cache>, kwargs...])

Emulate the continuous time evolution of given Hamiltonian `h` that has dependency on `t`. The default algorithm
we use is `Tsit5`. Available `kwargs` can be found in [Common Solver Options](https://diffeq.sciml.ai/latest/basics/common_solver_opts/).
Available algorithms can be found at [Full List of Methods](https://diffeq.sciml.ai/latest/solvers/ode_solve/#Full-List-of-Methods).

# Emulation Cache

Emulation of time evolution requires allocating intermediate variables that can be large. Thus preallocating
these intermediate variables and share this memory between iterations can speed up the simulation when `emulate!`
is called for multiple times for similar Hamiltonian. By similar Hamiltonian, we mean Hamiltonians that contains
the same term but can have different parameters.

## Time Independent Cache

The cache for time independent simulation contains two part: a) a `SparseMatrixCSC` matrix to store the intermediate
Hamiltonian matrices. b) a `KrylovSubspace` object to store the Krylov subspace.

## continuous Time Cache

For continuous time emulation, one only needs to create a `SparseMatrixCSC` matrix to store the intermediate
Hamiltonian matrices.

## Examples

The simplest way to use cache inside a large loop is via closure. The following example
returns a closure that takes time `ts` as parameters during simulation.

```jl
function your_emulation_task(n, subspace, hs)
    r = zero_state(n, subspace)
    cache = EmulatorCache(first(hs), subspace)
    return function task(ts)
        set_zero_state!(r)
        emulate!(r, ts, hs; cache=cache)
    end
end
```

or if you just want to use `ts` and `ϕs` as your parameters

```jl
function your_emulation_task(n, subspace)
    r = zero_state(n, subspace)
    cache = EmulatorCache(simple_rydberg(n, 1.0), subspace)

    return function task(xs)
        set_zero_state!(r)
        ts = xs[1:length(xs)÷2]
        ϕs = xs[length(xs)÷2+1:end]
        hs = simple_rydberg.(n, ϕs)
        emulate!(r, ts, hs; cache=cache)
    end
end
```
"""
function emulate! end

emulate_cache(r::Yao.ArrayReg, h::AbstractTerm) = EmulatorCache(real(eltype(r.state)), h)
emulate_cache(r::RydbergReg, h::AbstractTerm) = EmulatorCache(real(eltype(r.state)), h, r.subspace)

function emulate!(r::Yao.AbstractRegister, t::Number, h::AbstractTerm; cache=emulate_cache(r, h))
    emulate_routine!(r, default_unit(μs, t), h, cache.H)
    return r
end

function emulate!(r::Yao.AbstractRegister, ts::Vector{<:Number}, hs::Vector{<:AbstractTerm}; cache=emulate_cache(r, first(hs)))
    @progress for (t, h) in zip(ts, hs)
        emulate!(r, t, h; cache)
    end
    return r
end

"""
    emulate(t_or_ts, h_or_hs; kwargs...)

Non in-place version of [`emulate!`](@ref). See [`emulate!`](@ref) for valid kwargs.
"""
function emulate(t_or_ts, h_or_hs; kwargs...)
    precision_t = real(eltype(t_or_ts))
    r = Yao.zero_state(Complex{precision_t}, nsites(h_or_hs))
    return emulate!(r, t_or_ts, h_or_hs; kwargs...)
end

"""
    emulate(s::Subspace, t_or_ts, h_or_hs; kwargs...)

Non in-place version of [`emulate!`](@ref). See [`emulate!`](@ref) for valid kwargs.
"""
function emulate(s::Subspace, t_or_ts, h_or_hs; kwargs...)
    precision_t = real(eltype(t_or_ts))
    r = zero_state(Complex{precision_t}, nsites(h_or_hs), s)
    return emulate!(r, t_or_ts, h_or_hs; kwargs...)
end
