function qaoa_routine!(r::RydbergReg, ts::Vector, hs::Vector{<: AbstractTerm}, Ks::KrylovSubspace, cache)
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

function qaoa_routine!(r::Yao.ArrayReg, ts::Vector, hs::Vector{<: AbstractTerm}, Ks::KrylovSubspace, cache)
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

struct QAOA{K <: KrylovSubspace, C}
    Ks::K
    cache::C
end

function QAOA(N::Int, H::AbstractMatrix{Complex{T}}; krylov_niteration=min(30, N)) where T
    Ks = KrylovSubspace{Complex{T}, T}(N, krylov_niteration)
    return QAOA(Ks, H)
end

# TODO: verify the nonzero structure of given Hamiltonian sequence
function QAOA(::Type{T}, H::AbstractTerm, n::Int; kwargs...) where {T <: Real}
    return QAOA(1 << n, SparseMatrixCSC{Complex{T}}(H); kwargs...)
end

function QAOA(::Type{T}, H::AbstractTerm, s::Subspace; kwargs...) where {T <: Real}
    return QAOA(length(s), SparseMatrixCSC{Complex{T}}(H, s); kwargs...)
end

QAOA(H::AbstractTerm, xs...; kwargs...) = QAOA(Float64, H, xs...; kwargs...)

# TODO: make full use of YaoLang interface
function (u::QAOA)(ts::Vector, hs::Vector{<:AbstractTerm})
    function routine(r)
        qaoa_routine!(r, ts, hs, u.Ks, u.cache)
        return r
    end
end
