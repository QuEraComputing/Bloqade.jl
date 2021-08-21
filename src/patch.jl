# CUDA patch
using LinearAlgebra

# https://github.com/JuliaGPU/CUDA.jl/pull/217
# symmetric mul!
# level 2
@inline function LinearAlgebra.mul!(y::CuVector{T}, A::Hermitian{T,<:CuMatrix}, x::CuVector{T},
    α::Number, β::Number) where {T<:CUDA.CUBLAS.CublasReal}
    alpha, beta = promote(α, β, zero(T))
    if alpha isa Union{Bool,T} && beta isa Union{Bool,T}
        return CUBLAS.symv!(A.uplo, alpha, A.data, x, beta, y)
    else
        error("only supports BLAS type, got $T")
    end
end

@inline function LinearAlgebra.mul!(y::CuVector{T}, A::Hermitian{T,<:CuMatrix}, x::CuVector{T},
    α::Number, β::Number) where {T<:CUDA.CUBLAS.CublasComplex}
    alpha, beta = promote(α, β, zero(T))
    if alpha isa Union{Bool,T} && beta isa Union{Bool,T}
        return CUBLAS.hemv!(A.uplo, alpha, A.data, x, beta, y)
    else
        error("only supports BLAS type, got $T")
    end
end

# level 3
@inline function LinearAlgebra.mul!(C::CuMatrix{T}, A::Hermitian{T,<:CuMatrix}, B::CuMatrix{T},
    α::Number, β::Number) where {T<:CUDA.CUBLAS.CublasReal}
    alpha, beta = promote(α, β, zero(T))
    if alpha isa Union{Bool,T} && beta isa Union{Bool,T}
        return CUBLAS.symm!('L', A.uplo, alpha, A.data, B, beta, C)
    else
        error("only supports BLAS type, got $T")
    end
end

@inline function LinearAlgebra.mul!(C::CuMatrix{T}, A::CuMatrix{T}, B::Hermitian{T,<:CuMatrix},
    α::Number, β::Number) where {T<:CUDA.CUBLAS.CublasReal}
    alpha, beta = promote(α, β, zero(T))
    if alpha isa Union{Bool,T} && beta isa Union{Bool,T}
        return CUBLAS.symm!('R', B.uplo, alpha, B.data, A, beta, C)
    else
        error("only supports BLAS type, got $T")
    end
end

@inline function LinearAlgebra.mul!(C::CuMatrix{T}, A::Hermitian{T,<:CuMatrix}, B::CuMatrix{T},
    α::Number, β::Number) where {T<:CUDA.CUBLAS.CublasComplex}
    alpha, beta = promote(α, β, zero(T))
    if alpha isa Union{Bool,T} && beta isa Union{Bool,T}
        return CUBLAS.hemm!('L', A.uplo, alpha, A.data, B, beta, C)
    else
        error("only supports BLAS type, got $T")
    end
end

@inline function LinearAlgebra.mul!(C::CuMatrix{T}, A::CuMatrix{T}, B::Hermitian{T,<:CuMatrix},
    α::Number, β::Number) where {T<:CUDA.CUBLAS.CublasComplex}
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

# NOTE: this doesn't have much improvement, so we don't support it for now
# we also need to fix the problem of precision in upstream on our Hamiltonian matrix

# # CUDA.expv!
# function ExponentialUtilities.expv!(w::CuVector{Complex{Tw}}, t::Complex, Ks::KrylovSubspace{T, U};
#       cache=nothing, dexpHe::CuVector = CuVector{Complex{real(U)}}(undef, Ks.m)) where {Tw, T, U}
#     m, beta, V, H = Ks.m, Ks.beta, getV(Ks), getH(Ks)
#     @assert length(w) == size(V, 1) "Dimension mismatch"
#     if cache === nothing
#         cache = Matrix{U}(undef, m, m)
#     elseif isa(cache, ExpvCache)
#         cache = get_cache(cache, m)
#     else
#         throw(ArgumentError("Cache must be an ExpvCache"))
#     end
#     copyto!(cache, @view(H[1:m, :]))
#     if ishermitian(cache)
#         # Optimize the case for symtridiagonal H
#         F = eigen!(SymTridiagonal(cache))
#         expHe = F.vectors * (exp.(lmul!(t,F.values)) .* @view(F.vectors[1, :]))
#     else
#         expH = t .* cache # Complex
#         _exp!(expH)
#         expHe = @view(expH[:, 1])
#     end
#     copyto!(dexpHe, expHe)
#     lmul!(beta, mul!(w, @view(V[:, 1:m]), @view(dexpHe[1:m]))) # exp(A) ≈ norm(b) * V * exp(H)e
# end

# using Cassette, KernelAbstractions

# function Cassette.overdub(ctx::KernelAbstractions.CUDACtx, ::typeof(BitBasis.log2i), x::Int64)
#     @Base._inline_meta
#     return 63 - leading_zeros(x)
# end

# https://github.com/JuliaGPU/CUDA.jl/pull/1106

struct CuSparseDeviceMatrixCSR{Tv} <: AbstractCuSparseMatrix{Tv}
    rowPtr::CuDeviceVector{Cint, AS.Global}
    colVal::CuDeviceVector{Cint, AS.Global}
    nzVal::CuDeviceVector{Tv, AS.Global}
    dims::NTuple{2, Int}
    nnz::Cint
end

Base.size(H::CuSparseDeviceMatrixCSR) = H.dims

function Adapt.adapt_structure(to::CUDA.Adaptor, x::CuSparseMatrixCSR{Tv}) where Tv
    CuSparseDeviceMatrixCSR(
        cudaconvert(x.rowPtr),
        cudaconvert(x.colVal),
        cudaconvert(x.nzVal),
        x.dims, x.nnz
    )
end

function Base.show(io::IO, ::MIME"text/plain", A::CuSparseDeviceMatrixCSR)
    println(io, "$(length(A))-element device sparse matrix CSR at:")
    println(io, "  rowPtr $(pointer(A.rowPtr))")
    println(io, "  colVal $(pointer(A.colVal))")
    print(io, "  nzVal $(pointer(A.nzVal))")
end

function Adapt.adapt_structure(to, t::XTerm)
    XTerm(t.nsites, adapt(to, t.Ωs), adapt(to, t.ϕs))
end

function Adapt.adapt_structure(to, t::ZTerm)
    ZTerm(t.nsites, adapt(to, t.Δs))
end

function Adapt.adapt_structure(to, t::RydInteract)
    RydInteract(adapt(to, t.atoms), t.C)
end

function Adapt.adapt_structure(to, t::Hamiltonian)
    Hamiltonian(map(x->Adapt.adapt(to, x), t.terms))
end

function Adapt.adapt_structure(to, r::RydbergReg{N}) where {N}
    return RydbergReg{N}(adapt(to, r.state), adapt(to, r.subspace))
end

Adapt.adapt_structure(to, s::Subspace) = Subspace(s.map, adapt(to, s.subspace_v))

function Adapt.adapt_structure(to, cache::EmulatorCache)
    return EmulatorCache(CuSparseMatrixCSR(cache.H))
end
