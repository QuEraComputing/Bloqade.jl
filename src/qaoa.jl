using Yao
using ExponentialUtilities
using SparseArrays

# NOTE: we use a type parameter Hs here to let Julia specialize the type when possible.
# but this should be a Vector in general, since we need a mutable array object
struct QAOA{N, Hs <: Vector, TimeType <: Real, HMatrixType <: AbstractMatrix{Complex{TimeType}},KrylovT <: KrylovSubspace} <: Yao.AbstractBlock{N}
    hamiltonians::Hs
    hamiltonian_cache::HMatrixType
    ts::Vector{TimeType}
    # we use Int for now since we are not
    # targeting system above 64 qubits for now
    subspace_v::Vector{Int}
    Ks::KrylovT
end

function QAOA{N}(subspace_v::Vector{Int}, hs::Vector, ts::Vector{TimeType};
        cache=(m=length(subspace_v);spzeros(Complex{TimeType}, m, m)),
        krylov_m=min(30, length(subspace_v))) where {N, TimeType}
    
    m = length(subspace_v)
    Ks = KrylovSubspace{Complex{TimeType}, TimeType}(m, krylov_m)
    QAOA{N, typeof(hs), typeof(ts), typeof(cache), typeof(Ks)}(hs, cache, ts, subspace_v, Ks)
end

# only defined for single batch for now
function Yao.apply!(r::AbstractRegister{1}, x::QAOA{N, Hs, T}) where {N, Hs, T}
    for (h, t) in zip(x.hamiltonians, x.ts)
        to_matrix!(x.hamiltonian_cache, N, x.subspace_v, one(T), h.ϕ)
        # qaoa step
        # NOTE: we share the Krylov subspace here since
        #       the Hamiltonians have the same shape
        arnoldi!(x.Ks, x.hamiltonian_cache, r.state; m=Ks_m, ishermitian=true)
        st = expv!(st, -im*t, Ks)
        fill!(H, zero(Complex{T}))
    end
    return st
end

"""
    evaluate_qaoa!(st::Vector{Complex{T}}, hs::Vector{<:AbstractRydbergHamiltonian}, n, subspace_v, ts::Vector{<:Real})

Evaluate a QAOA sequence `hs` along with parameters `ts` given initial state `st` and atom geometry `atoms`.
"""
function evaluate_qaoa! end

function evaluate_qaoa!(st::Vector{Complex{T}}, hs::Vector{SimpleRydberg{T}}, n::Int, subspace_v, ts::Vector{T}) where T
    m = length(subspace_v)
    H = spzeros(Complex{T}, m, m)

    # Krylov Subspace Cfg
    Ks_m = min(30, size(H, 1))
    Ks = KrylovSubspace{Complex{T}, T}(length(st), Ks_m)

    for (h, t) in zip(hs, ts)
        to_matrix!(H, n, subspace_v, one(T), h.ϕ)
        # qaoa step
        # NOTE: we share the Krylov subspace here since
        #       the Hamiltonians have the same shape
        arnoldi!(Ks, H, st; m=Ks_m, ishermitian=true)
        st = expv!(st, -im*t, Ks)
        dropzeros!(fill!(H, zero(Complex{T})))
    end
    return st
end
