using CUDA
using ExponentialUtilities
using LinearAlgebra
using ExponentialUtilities: getV, getH, get_cache, _exp!
using LinearAlgebra: BlasReal, BlasComplex
using SparseArrays
using CUDA: CUBLAS
using CUDA: GPUArrays
using CUDA.GPUArrays: AbstractGPUVecOrMat, AbstractGPUArray, AbstractGPUVector

function CUDA.cu(r::RydbergReg{N}) where {N}
    return RydbergReg{N}(cu(r.state), cu(r.subspace))
end

# CUDA.allowscalar(false)

# CUDA patch
# symmetric mul!
# level 2
@inline function LinearAlgebra.mul!(y::CuVector{T}, A::Hermitian{T,<:CuMatrix}, x::CuVector{T},
             α::Number, β::Number) where {T<:BlasReal}
    alpha, beta = promote(α, β, zero(T))
    if alpha isa Union{Bool,T} && beta isa Union{Bool,T}
        return CUBLAS.symv!(A.uplo, alpha, A.data, x, beta, y)
    else
        error("only supports BLAS type, got $T")
    end
end

@inline function LinearAlgebra.mul!(y::CuVector{T}, A::Hermitian{T,<:CuMatrix}, x::CuVector{T},
             α::Number, β::Number) where {T<:BlasComplex}
    alpha, beta = promote(α, β, zero(T))
    if alpha isa Union{Bool,T} && beta isa Union{Bool,T}
        return CUBLAS.hemv!(A.uplo, alpha, A.data, x, beta, y)
    else
        error("only supports BLAS type, got $T")
    end
end

# level 3

@inline function LinearAlgebra.mul!(C::CuMatrix{T}, A::Hermitian{T,<:CuMatrix}, B::CuMatrix{T},
             α::Number, β::Number) where {T<:BlasReal}
    alpha, beta = promote(α, β, zero(T))
    if alpha isa Union{Bool,T} && beta isa Union{Bool,T}
        return CUBLAS.symm!('L', A.uplo, alpha, A.data, B, beta, C)
    else
        error("only supports BLAS type, got $T")
    end
end
@inline function LinearAlgebra.mul!(C::CuMatrix{T}, A::CuMatrix{T}, B::Hermitian{T,<:CuMatrix},
             α::Number, β::Number) where {T<:BlasReal}
    alpha, beta = promote(α, β, zero(T))
    if alpha isa Union{Bool,T} && beta isa Union{Bool,T}
        return CUBLAS.symm!('R', B.uplo, alpha, B.data, A, beta, C)
    else
        error("only supports BLAS type, got $T")
    end
end
@inline function LinearAlgebra.mul!(C::CuMatrix{T}, A::Hermitian{T,<:CuMatrix}, B::CuMatrix{T},
             α::Number, β::Number) where {T<:BlasComplex}
    alpha, beta = promote(α, β, zero(T))
    if alpha isa Union{Bool,T} && beta isa Union{Bool,T}
        return CUBLAS.hemm!('L', A.uplo, alpha, A.data, B, beta, C)
    else
        error("only supports BLAS type, got $T")
    end
end
@inline function LinearAlgebra.mul!(C::CuMatrix{T}, A::CuMatrix{T}, B::Hermitian{T,<:CuMatrix},
             α::Number, β::Number) where {T<:BlasComplex}
    alpha, beta = promote(α, β, zero(T))
    if alpha isa Union{Bool,T} && beta isa Union{Bool,T}
        return CUBLAS.hemm!('R', B.uplo, alpha, B.data, A, beta, C)
    else
        error("only supports BLAS type, got $T")
    end
end

# GPUArray patch
@inline function LinearAlgebra.mul!(C::AbstractGPUArray{Complex{Float64},1},
        A::AbstractGPUVecOrMat{Complex{Float64}}, B::AbstractGPUVector, a::Real, b::Real)
    GPUArrays.generic_matmatmul!(C, A, B, a, b)
end

# CUDA.expv!

function ExponentialUtilities.expv!(w::CuVector{Tw}, t::Real, Ks::KrylovSubspace{T, U};
               cache=nothing, dexpHe::CuVector = CuVector{U}(undef, Ks.m)) where {Tw, T, U}
    m, beta, V, H = Ks.m, Ks.beta, getV(Ks), getH(Ks)
    @assert length(w) == size(V, 1) "Dimension mismatch"
    if cache === nothing
        cache = Matrix{U}(undef, m, m)
    elseif isa(cache, ExpvCache)
        cache = get_cache(cache, m)
    else
        throw(ArgumentError("Cache must be an ExpvCache"))
    end
    copyto!(cache, @view(H[1:m, :]))
    if ishermitian(cache)
        # Optimize the case for symtridiagonal H
        F = eigen!(SymTridiagonal(cache))
        expHe = F.vectors * (exp.(lmul!(t,F.values)) .* @view(F.vectors[1, :]))
    else
        lmul!(t, cache); expH = cache
        _exp!(expH)
        expHe = @view(expH[:, 1])
    end

    copyto!(dexpHe, expHe)
    lmul!(beta, mul!(w, @view(V[:, 1:m]), @view(dexpHe[1:m]))) # exp(A) ≈ norm(b) * V * exp(H)e
end


function ExponentialUtilities.expv!(w::CuVector{Complex{Tw}}, t::Real, Ks::KrylovSubspace{T, U};
        cache=nothing, dexpHe::CuVector = CuVector{U}(undef, Ks.m)) where {Tw, Tt, T, U}
    m, beta, V, H = Ks.m, Ks.beta, getV(Ks), getH(Ks)
    @assert length(w) == size(V, 1) "Dimension mismatch"
    if cache === nothing
        cache = Matrix{U}(undef, m, m)
    elseif isa(cache, ExpvCache)
        cache = get_cache(cache, m)
    else
        throw(ArgumentError("Cache must be an ExpvCache"))
    end
    copyto!(cache, @view(H[1:m, :]))
    if ishermitian(cache)
        # Optimize the case for symtridiagonal H
        F = eigen!(SymTridiagonal(cache))
        expHe = F.vectors * (exp.(t * F.values) .* @view(F.vectors[1, :]))
    else
        lmul!(t, cache); expH = cache
        _exp!(expH)
        expHe = @view(expH[:, 1])
    end

    copyto!(dexpHe, expHe)
    lmul!(beta, mul!(w, @view(V[:, 1:m]), @view(dexpHe[1:m]))) # exp(A) ≈ norm(b) * V * exp(H)e
end

export CuQAOA
# hamiltonian
struct CuQAOA{N, Hs <: Vector, TimeType <: Real, KrylovT <: KrylovSubspace} <: Yao.PrimitiveBlock{N}
    hamiltonians::Hs
    ts::Vector{TimeType}
    # we use Int for now since we are not
    # targeting system above 64 qubits for now
    subspace_v::Vector{Int}
    device_subspace_v::CuVector{Int}
    Ks::KrylovT

    # cache
    # TODO: make a cache type and unify QAOA struct
    h_cache::CuSparseMatrixCSR{Complex{TimeType}}
    dexpHe::CuVector{TimeType}
end

function CuQAOA{N}(subspace_v::Vector{Int}, hs::Vector, ts::Vector{TimeType};
        cache=nothing,
        krylov_niteration=min(30, length(subspace_v))) where {N, TimeType}

    m = length(subspace_v)
    V = CuMatrix{Complex{TimeType}}(undef, m, krylov_niteration + 1)
    H = fill(zero(TimeType), krylov_niteration + 1, krylov_niteration)
    Ks = KrylovSubspace{Complex{TimeType}, TimeType, TimeType, typeof(V), typeof(H)}(krylov_niteration, krylov_niteration, false, zero(TimeType), V, H);
    if cache === nothing
        cpu_cache = init_hamiltonian(TimeType, N, m, subspace_v, one(TimeType), first(hs).ϕ)
        m_cache = CuSparseMatrixCSR(cpu_cache)
        dexpHe = CuVector{TimeType}(undef, Ks.maxiter)
    else
        m_cache, dexpHe = cache
    end
    return CuQAOA{N, typeof(hs), eltype(ts), typeof(Ks)}(hs, ts, subspace_v, cu(subspace_v), Ks, m_cache, dexpHe)
end

function Yao.apply!(r::RydbergReg{N, 1, <:CuArray, <:CuArray}, x::CuQAOA{N, Hs, T}) where {N, Hs, T}
    st = vec(r.state)
    qaoa_routine!(vec(r.state), x.hamiltonians, N, x.device_subspace_v, x.ts, x.Ks, x.h_cache, x.dexpHe)
    return r
end

# TODO: polish this routine so we don't need to have two similar routines
function qaoa_routine!(st::CuVector{Complex{T}}, hs::Vector{<:SimpleRydberg}, n::Int,
        subspace_v::CuVector, ts::Vector{T}, Ks::KrylovSubspace, cache::CuSparseMatrixCSR, dexpHe::CuVector) where T

    for (h, t) in zip(hs, ts)
        update_hamiltonian!(cache, n, subspace_v, one(T), h.ϕ)
        # qaoa step
        # NOTE: we share the Krylov subspace here since
        #       the Hamiltonians have the same shape
        lmul!(-im, cache.nzVal)
        arnoldi!(Ks, cache, st; ishermitian=false)
        expv!(st, t, Ks; dexpHe=dexpHe)
    end
    return st
end

function update_hamiltonian!(dst::CuSparseMatrixCSR, n::Int, subspace_v::CuVector, Ω, ϕ)
    function kernel(rowPtr, colVal, nzVal, subspace_v, Ω, ϕ)
        row = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        if row > length(rowPtr) - 1
            return
        end

        for k in rowPtr[row]:rowPtr[row+1]-1
            col = colVal[k]
            lhs = subspace_v[row]
            rhs = subspace_v[col]

            mask = lhs ⊻ rhs
            if lhs & mask == 0
                nzVal[k] = Ω * exp(im * ϕ)
            else
                nzVal[k] = Ω * exp(-im * ϕ)
            end
            # update_x_term!(nzVal, lhs, rhs, k, Ω, ϕ)
        end
        return
    end
    
    if size(dst, 2) < 256
        threads = size(dst, 2)
        nblocks = 1
    else
        threads = 256
        nblocks = ceil(Int, size(dst, 2) / 256)
    end

    @cuda threads=threads blocks=nblocks kernel(dst.rowPtr, dst.colVal, dst.nzVal, subspace_v, Ω, ϕ)
    return dst
end

function update_hamiltonian!(dst::CuSparseMatrixCSR, n::Int, subspace_v::CuVector, Ω, ϕ, Δ)
    function kernel(rowPtr, colVal, nzVal, subspace_v, Ω, ϕ, Δ)
        row = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        if row > length(rowPtr) - 1
            return
        end

        for k in rowPtr[row]:rowPtr[row+1]-1
            col = colVal[k]
            lhs = subspace_v[row]
            rhs = subspace_v[col]

            if row == col
                update_z_term!(nzVal, lhs, n, k, Δ)
            else
                update_x_term!(nzVal, lhs, rhs, k, Ω, ϕ)
            end
        end
        return
    end

    if size(dst, 2) < 256
        threads = size(dst, 2)
        nblocks = 1
    else
        threads = 256
        nblocks = ceil(Int, size(dst, 2) / 256)
    end

    @cuda threads=threads blocks=nblocks kernel(dst.rowPtr, dst.colVal, dst.nzVal, subspace_v, Ω, ϕ, Δ)
    return dst
end
