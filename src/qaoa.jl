using Yao
using ExponentialUtilities
using SparseArrays
using LuxurySparse

export QAOA, update_ansatz!
# NOTE: we use a type parameter Hs here to let Julia specialize the type when possible.
# but this should be a Vector in general, since we need a mutable array object
struct QAOA{N, Hs <: Vector, TimeType <: Real, HMatrixType <: AbstractMatrix{Complex{TimeType}},KrylovT <: KrylovSubspace} <: Yao.PrimitiveBlock{N}
    hamiltonians::Hs
    hamiltonian_cache::HMatrixType
    ts::Vector{TimeType}
    # we use Int for now since we are not
    # targeting system above 64 qubits for now
    subspace_v::Vector{Int}
    Ks::KrylovT
end

function init_hamiltonian(::Type{T}, n::Int, m::Int, subspace_v, Ω, ϕ) where T
    H = SparseMatrixCOO{Complex{T}}(undef, m, m)
    to_matrix!(H, n, subspace_v, Ω, ϕ)
    return SparseMatrixCSC(H)
end

"""
    QAOA{N}(subspace_v, hs::Vector, ts::Vector{TimeType}; kwargs...)

Create a QAOA block, where `subspace_v` is the emulation subspace, `hs` is the vector of hamiltonians,
`ts` is the time parameters.

### Keywords

- `cache=(m=length(subspace_v);spzeros(Complex{TimeType}, m, m))`, hamiltonian matrix cache
- `krylov_niteration=min(30, length(subspace_v))`, maximum number of iteration for Krylov method
"""
function QAOA{N}(subspace_v::Vector{Int}, hs::Vector, ts::Vector{TimeType};
        cache=nothing,
        krylov_niteration=min(30, length(subspace_v))) where {N, TimeType}
    
    m = length(subspace_v)
    Ks = KrylovSubspace{Complex{TimeType}, TimeType}(m, krylov_niteration)

    if cache === nothing
        cache = init_hamiltonian(TimeType, N, m, subspace_v, one(TimeType), first(hs).ϕ)
    end
    return QAOA{N, typeof(hs), eltype(ts), typeof(cache), typeof(Ks)}(hs, cache, ts, subspace_v, Ks)
end

# TODO: move this to utils.jl
"""
    fillzero!(M)

Fill a given array `M` with zero inplace.
"""
fillzero!(M::AbstractArray) = fill!(M, 0)

function fillzero!(M::SparseMatrixCSC)
    fill!(M.colptr, 1)
    resize!(M.nzval, 0)
    resize!(M.rowval, 0)
    return M
end

# only defined for single batch for now
function Yao.apply!(r::RydbergReg{N, 1}, x::QAOA{N, Hs, T}) where {N, Hs, T}
    st = vec(r.state)
    qaoa_routine!(vec(r.state), x.hamiltonians, N, x.subspace_v, x.ts, x.Ks, x.hamiltonian_cache)
    return r
end

# NOTE: very basic interface to update SimpleRydberg QAOA for now
# need to specialize for other kinds of hamiltonians or make it
# more generic in the future
#
# NOTE: RL: I use update_ansatz! instead of update! here to distinguish with Flux.Optimise.update!
#           so we can use Flux.Optimise for Policy Gradients.
function update_ansatz!(x::QAOA{N, Vector{SimpleRydberg{T}}, T}, ϕ::Vector{T}, ts::Vector{T}) where {N, T}
    @assert length(ϕ) == length(ts) == length(x.hamiltonians) "number of input parameters does not match ansatz parameters"
    @inbounds for k in 1:length(ϕ)
        x.hamiltonians[k] = SimpleRydberg(ϕ[k])
        x.ts[k] = ts[k]
    end
    return x
end

function qaoa_routine!(st::Vector{Complex{T}}, hs::Vector{SimpleRydberg{T}}, n::Int, subspace_v, ts::Vector{T}, Ks::KrylovSubspace, cache::AbstractMatrix) where T
    for (h, t) in zip(hs, ts)
        to_matrix!(cache, n, subspace_v, one(T), h.ϕ)
        # qaoa step
        # NOTE: we share the Krylov subspace here since
        #       the Hamiltonians have the same shape
        arnoldi!(Ks, cache, st)
        st = expv!(st, -im*t, Ks)
        fillzero!(cache)
    end
    return st
end
