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
    emulate!(r, ts, hs[, cache=EmulatorCache(ts, hs)])

Emulate the time evolution of a sequence of Hamiltonians `hs` of time length `ts` for each Hamiltonian.
An optional argument can be feeded to preallocate the memory.
"""
function emulate! end

function emulate!(r::Yao.ArrayReg, ts::Vector{<:Real}, hs::Vector{<:AbstractTerm}, cache=EmulatorCache(ts, hs))
    emulate_routine!(r, ts, hs, cache.Ks, cache.H)
    return r
end

function emulate!(r::RydbergReg, ts::Vector{<:Real}, hs::Vector{<:AbstractTerm}, cache=EmulatorCache(ts, hs, r.subspace))
    emulate_routine!(r, ts, hs, cache.Ks, cache.H)
    return r
end
