function CUDA.cu(r::RydbergReg{N}) where {N}
    return RydbergReg{N}(cu(r.state), cu(r.subspace))
end

using CUDA
using ExponentialUtilities
using LinearAlgebra
using ExponentialUtilities: getV, getH, get_cache, _exp!
using LinearAlgebra: BlasReal, BlasComplex
using SparseArrays
using CUDA: CUBLAS
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

# CUDA.expv!

function ExponentialUtilities.expv!(w::CuVector{Tw}, t::Real, Ks::KrylovSubspace{T, U};
               cache=nothing, dexpHe::CuVector = CuVector{U}(undef, Ks.m)) where {Tw, T, U}
    m, beta, V, H = Ks.m, Ks.beta, getV(Ks), getH(Ks)
    @assert length(w) == size(V, 1) "Dimension mismatch"
    if cache == nothing
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
    lmul!(beta, mul!(w, @view(V[:, 1:m]), dexpHe)) # exp(A) ≈ norm(b) * V * exp(H)e
end

# hamiltonian
function update_z_term!(dst::CuSparseMatrixCSR{T}, n::Int, count::Int, col::Int, lhs::Int, Δ) where T
    sigma_z = zero(T)
    for k in 1:n
        if readbit(lhs, k) == 1
            sigma_z -= getscalarmaybe(Δ, k)
        else
            sigma_z += getscalarmaybe(Δ, k)
        end
    end
    dst.nzval[count] = sigma_z
    return dst
end

function update_x_term!(dst::SparseMatrixCSC, count::Int, row::Int, col::Int, lhs::Int, rhs::Int, Ω, ϕ)
    mask = lhs ⊻ rhs
    k = log2i(mask) + 1
    if lhs & mask == 0
        dst.nzval[count] = getscalarmaybe(Ω, k) * exp(im * getscalarmaybe(ϕ, k))
    else
        dst.nzval[count] = getscalarmaybe(Ω, k) * exp(-im * getscalarmaybe(ϕ, k))
    end
    return dst
end

function update_hamiltonian!(dst::SparseMatrixCSC, n::Int, subspace_v, Ω, ϕ)
    col = 1
    for (count, v) in enumerate(dst.nzval)
        if count == dst.colptr[col+1]
            col += 1
        end

        row = dst.rowval[count]
        lhs = subspace_v[row]
        rhs = subspace_v[col]
        # we don't check if row == col
        # since there is only x term
        # update x term
        update_x_term!(dst, count, row, col, lhs, rhs, Ω, ϕ)
    end
    return dst
end
