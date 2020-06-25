struct QAOA{N, T, Term, C, S, K, E} <: Yao.AbstractBlock{N, T}
    ts::Vector{T}
    term::Vector{Term}
    cache::C
    subspace::S
    Ks::K
    expHe::E # expv cache
end

"""
    QAOA{N}(subspace_v, hs::Vector, ts::Vector{TimeType}; kwargs...)

Create a QAOA block, where `subspace_v` is the emulation subspace, `hs` is the vector of hamiltonians,
`ts` is the time parameters.

### Keywords

- `cache=(m=length(subspace_v);spzeros(Complex{TimeType}, m, m))`, hamiltonian matrix cache
- `krylov_niteration=min(30, length(subspace_v))`, maximum number of iteration for Krylov method
"""
function QAOA{N}(ts::Vector{T}, hs::Vector{Term};
        cache=nothing, subspace=nothing, krylov_niteration=min(30, length(subspace))) where {N, T, Term}

    if subspace === nothing
        m = nsites(hs[1])
    else
        m = length(subspace)
    end

    Ks = KrylovSubspace{Complex{T}, T}(m, krylov_niteration)

    if cache === nothing
        H = SparseMatrixCOO{Complex{T}}(undef, m, m)
        if subspace === nothing
            to_matrix!(H, hs[1])
        else
            to_matrix!(H, hs[1], subspace)
        end
        cache = SparseMatrixCSC(H)
    end

    C = typeof(cache)
    S = typeof(subspace)
    return QAOA{N, T, Term, C, S, typeof(Ks), Nothing}(ts, hs, cache, subspace, Ks, nothing)
end

function Yao.apply!(r::RydbergReg{N, 1}, x::QAOA{N, T}) where {N, T}
    st = vec(r.state)
    for (h, t) in zip(x.term, x.ts)
        update_term!(x.cache, h, x.subspace)
        # qaoa step
        # NOTE: we share the Krylov subspace here since
        #       the Hamiltonians have the same shape
        arnoldi!(x.Ks, x.cache, st)
        expv!(st, -im*t, x.Ks)
    end
    return r
end
