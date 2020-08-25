function emulate_routine!(r::RydbergReg, ts::Vector, hs::Vector{<: AbstractTerm}, Ks::KrylovSubspace, cache)
    st = vec(r.state)
    for (h, t) in zip(hs, ts)
        update_term!(cache, h, r.subspace)
        # qaoa step
        # NOTE: we share the Krylov subspace here since
        #       the Hamiltonians have the same shape
        arnoldi!(Ks, cache, st)
        expv!(st, -im*t, Ks)
    end
    return r
end

function emulate_routine!(r::Yao.ArrayReg, ts::Vector, hs::Vector{<: AbstractTerm}, Ks::KrylovSubspace, cache)
    st = vec(r.state)
    for (h, t) in zip(hs, ts)
        update_term!(cache, h)
        # qaoa step
        # NOTE: we share the Krylov subspace here since
        #       the Hamiltonians have the same shape
        arnoldi!(Ks, cache, st)
        expv!(st, -im*t, Ks)
    end
    return r
end

"""
    EmulatorCache

Cache type for the emulator.
"""
struct EmulatorCache{K <: KrylovSubspace, C}
    Ks::K
    H::C
end

function EmulatorCache(N::Int, H::AbstractMatrix{Complex{T}}; krylov_niteration=min(30, N)) where T
    Ks = KrylovSubspace{Complex{T}, T}(N, krylov_niteration)
    return EmulatorCache(Ks, H)
end

function EmulatorCache(::Type{T}, H::AbstractTerm, n::Int; kwargs...) where {T <: Real}
    return EmulatorCache(1 << n, SparseMatrixCSC{Complex{T}}(H); kwargs...)
end

function EmulatorCache(::Type{T}, H::AbstractTerm, s::Subspace; kwargs...) where {T <: Real}
    return EmulatorCache(length(s), SparseMatrixCSC{Complex{T}}(H, s); kwargs...)
end

EmulatorCache(H::AbstractTerm, xs...; kwargs...) = EmulatorCache(Float64, H, xs...; kwargs...)

"""
    emulate!(r, ts, hs[; cache=EmulatorCache(ts, hs)])

Emulate the time evolution of a sequence of Hamiltonians `hs` of time length `ts` for each Hamiltonian.
An optional argument can be feeded to preallocate the memory.

    emulate!(r, t::Real, h::AbstractTerm[; algo=Tsit5(), cache=<default cache>, kwargs...])

Emulate the contiguous time evolution of given Hamiltonian `h` that has dependency on `t`. The default algorithm
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

## Continouns Time Cache

For continouns time emulation, one only needs to create a `SparseMatrixCSC` matrix to store the intermediate
Hamiltonian matrices.
"""
function emulate! end

function emulate!(r::Yao.ArrayReg, ts::Vector{T}, hs::Vector{<:AbstractTerm}; cache=EmulatorCache(T, first(hs), nsites(first(hs)))) where {T <: Real}
    emulate_routine!(r, ts, hs, cache.Ks, cache.H)
    return r
end

function emulate!(r::RydbergReg, ts::Vector{T}, hs::Vector{<:AbstractTerm}; cache=EmulatorCache(T, first(hs), r.subspace)) where {T <: Real}
    emulate_routine!(r, ts, hs, cache.Ks, cache.H)
    return r
end
