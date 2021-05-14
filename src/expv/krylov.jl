mutable struct KrylovSubspace{VEType, HEType, T <: Real, VType <: AbstractMatrix{VEType}, HType <: AbstractMatrix{HEType}, View, HmView, Algo}
    krylov_dim::Int # subspace size
    tol::T # tolerance
    V::VType # orthonormal bases
    H::HType # Gram-Schmidt coefficients + 2 error estimation dim
    v::View # pre-created views of V
    Hm::HmView # H[1:m, 1:m]
    algo::Algo
end

function KrylovSubspace(t, A, b, krylov_dim::Int=min(30, size(A, 1)), tol::Real=1e-7)
    T = promote_expv(t, A, b)
    n = length(b)
    HEType = real(T)
    V = similar(b, T, (n, krylov_dim+1))
    views = [view(V, :, k) for k in 1:krylov_dim+1]
    
    # NOTE:
    # this matrix should always be a Matrix
    # for CPU calculation since the arnoldi/lanczos
    # process cannot be implemented efficiently
    # on GPU
    # algo = ishermitian(A) ? Lanczos() : Arnoldi()
    algo = Arnoldi()
    H = zeros(T, (krylov_dim+2, krylov_dim+2))
    Hm = view(H, 1:krylov_dim, 1:krylov_dim)
    return KrylovSubspace(krylov_dim, tol, V, H, views, Hm, algo)
end

struct Arnoldi end
struct Lanczos end

function krylov_subspace!(Ks::KrylovSubspace, A, b::AbstractVector, beta::Real=norm(b))
    # initialize
    fill!(Ks.H, zero(eltype(Ks.H)))
    iszero(beta) && return Ks
    @. Ks.v[1] = b / beta
    krylov_subspace_iteration!(Ks, A, Ks.algo)
    # this extra dim is for error estimation
    Ks.H[Ks.krylov_dim+2, Ks.krylov_dim+1] = one(eltype(Ks.H))
    return Ks
end

function krylov_subspace_iteration!(Ks::KrylovSubspace, A, ::Arnoldi)
    # arnoldi steps
    @inbounds for j in 1:Ks.krylov_dim
        x, y = Ks.v[j], Ks.v[j+1]
        # Ks.V[j+1] = A * Ks.V[j]
        mul!(y, A, x)
        for i in 1:j
            α = dot(Ks.v[i], y)
            # H[i, j] = V[i]' ⋅ V[j+1]
            Ks.H[i, j] = α
            # V[j+1] = V[j+1] - H[i, j] V[i]
            axpy!(-α, Ks.v[i], y)
        end
        β = norm(y)

        # happy breakdown
        if β < Ks.tol
            Ks.krylov_dim = j
            break
        end

        Ks.H[j+1, j] = β
        lmul!(1/β, y)
    end
    return Ks
end

# TODO: fix this
# function krylov_subspace_iteration!(Ks::KrylovSubspace, A, ::Lanczos)
#     # initial iteration step
#     # v_2 = A ⋅ v_1
#     # α = dot(w, v_1)
#     # v_2 = v_2 - α ⋅ v_1
#     v_1 = Ks.v[1]
#     v_2 = mul!(Ks.v[2], A, v_1)
#     Ks.H[1, 1] = α = dot(v_2, v_1)
#     axpy!(-α, v_1, v_2)

#     @inbounds for k in 2:Ks.krylov_dim
#         # β_k = norm(v_k)
#         # if β_k != 0
#         #    v_k = v_k/β_k
#         # else
#         #   pick a random v with euclidean norm 1
#         # end
#         #
#         # v_{k+1} = A ⋅ v_k
#         # α[k] = dot(v_{k+1}, v_k)
#         # v_{k+1} = v_{k+1} - α[k] v_k - β_k v_{k-1}
#         β = norm(Ks.v[k])
#         Ks.H[k-1, k] = Ks.H[k, k-1] = β
#         @show β
#         lmul!(1/β, Ks.v[k])

#         mul!(Ks.v[k+1], A, Ks.v[k])
#         α = dot(Ks.v[k+1], Ks.v[k])
#         Ks.H[k, k] = α

#         axpy!(-β, Ks.v[k-1], Ks.v[k+1])
#         axpy!(-α, Ks.v[k], Ks.v[k+1])

#         if Ks.tol > β
#             # happy-breakdown
#             Ks.krylov_dim = k
#             break
#         end
#     end
#     return Ks
# end
