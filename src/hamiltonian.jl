abstract type AbstractTerm end

struct RydInteract{T, Atom <: RydAtom} <: AbstractTerm
    C::T
    atoms::Vector{Atom}
end

struct XTerm{Omega, Phi} <: AbstractTerm
    Ωs::Omega
    ϕs::Phi
end

struct ZTerm{Delta} <: AbstractTerm
    Δs::Delta
end

struct Hamiltonian{Terms <: Tuple} <: AbstractTerm
    terms::Terms
end

XTerm(Ωs) = XTerm(Ωs, nothing)

"""
    nsites(term)

Return the number of sites of given Hamiltonian term.
"""
function nsites end

nsites(t::XTerm) = length(t.Ωs)
nsites(t::ZTerm) = length(t.Δs)
nsites(t::Hamiltonian) = nsites(t.terms[1])
nsites(t::RydInteract) = length(t.atoms)

hilbert_space(n::Int) = 0:((1<<n)-1)
hilbert_space(t::AbstractTerm) = hilbert_space(nsites(t))

function Base.show(io::IO, t::RydInteract)
    print(io, "C/|r_i - r_j|^6 ")
    printstyled(io, "n_i n_j", color=:light_blue)
end

function Base.show(io::IO, t::XTerm)
    print(io, "Ω ⋅ (e^{iϕ}")
    printstyled(io, "|0)⟨1|", color=:light_blue)
    print(io, "+e^{-iϕ}")
    printstyled(io, "|1⟩⟨0|", color=:light_blue)
end

function Base.show(io::IO, t::XTerm{<:Any, Nothing})
    print(io, "Ω")
    printstyled(io, " σ^x", color=:light_blue)
end

function Base.show(io::IO, t::ZTerm)
    print(io, "Δ")
    printstyled(io, " σ^z", color=:light_blue)
end

function Base.show(io::IO, x::Hamiltonian)
    print(io, x.terms[1])
    for t in x.terms[2:end]
        print(io, " + ")
        print(io, t)
    end
end

Base.:(+)(x::AbstractTerm, y::AbstractTerm) = Hamiltonian((x, y))
Base.:(+)(x::AbstractTerm, y::Hamiltonian) = Hamiltonian((x, y.terms...))
Base.:(+)(x::Hamiltonian, y::AbstractTerm) = Hamiltonian((x.terms..., y))
Base.:(+)(x::Hamiltonian, y::Hamiltonian) = Hamiltonian((x.terms..., y.terms...))

"""
    getterm(terms, k, k_site)

Get the value of k-th local term in `terms`
given the site configuration as `k_site`.
"""
function getterm end

function getterm(t::XTerm, k, k_site)
    if k_site == 0
        return getscalarmaybe(t.Ωs, k) * exp(im * getscalarmaybe(t.ϕs, k))
    else
        return getscalarmaybe(t.Ωs, k) * exp(-im * getscalarmaybe(t.ϕs, k))
    end
end

function getterm(t::ZTerm, k, k_site)
    if k_site == 0
        return getscalarmaybe(t.Δs, k)
    else
        return -getscalarmaybe(t.Δs, k)
    end
end

function getterm(t::Hamiltonian, k, k_site)
    error("composite Hamiltonian term cannot be indexed")
end

"""
    to_matrix!(dst, term[, subspace])

Create given term to a matrix and assign it to `dst`. An optional argument `subspace`
can be taken to construct the matrix in subspace. `dst` should be initialized to zero
entries by the user.
"""
function to_matrix! end

to_matrix(term::AbstractTerm, xs...) = to_matrix(ComplexF64, term, xs...)

function to_matrix(::Type{T}, term::AbstractTerm) where T
    N = 1 << nsites(term)
    return to_matrix!(spzeros(T, N, N), term)
end

function to_matrix(::Type{T}, term::AbstractTerm, s::Subspace) where T
    N = length(s)
    return to_matrix!(spzeros(T, N, N), term, s)
end

# full space
# C/|r_i - r_j|^6 n_i n_j
function to_matrix!(dst::AbstractMatrix{T}, t::RydInteract) where T
    for i in 1:nsites(t), j in 1:nsites(t)
        r_i, r_j = t.atoms[i], t.atoms[j]
        alpha = t.C / norm(r_i - r_j)^6
        # |11⟩⟨11|
        for lhs in itercontrol(nsites(t), [i, j], [1, 1])
            dst[lhs, lhs] = alpha
        end
    end
    return dst
end

# Ω ⋅ (e^{iϕ}|0)⟨1| + e^{-iϕ} |1⟩⟨0|)
function to_matrix!(dst::AbstractMatrix{T}, t::XTerm) where T
    @inbounds for lhs in hilbert_space(t)
        for k in 1:nsites(t)
            k_site = readbit(lhs, k)
            rhs = flip(lhs, 1 << (k - 1))
            dst[lhs+1, rhs+1] = getterm(t, k, k_site)
        end
    end
    return dst
end

function to_matrix!(dst::AbstractMatrix{T}, t::ZTerm) where T
    @inbounds for lhs in hilbert_space(t)
        sigma_z = zero(T)
        for k in 1:nsites(t)
            sigma_z += getterm(t, k, readbit(lhs, k))
        end
        dst[lhs+1, lhs+1] = sigma_z
    end
    return dst
end

function to_matrix!(dst::AbstractMatrix{T}, t::Hamiltonian, xs...) where T
    for term in t.terms
        to_matrix!(dst, term, xs...)
    end
    return dst
end

# subspace
function to_matrix!(dst::AbstractMatrix, t::XTerm, s::Subspace)
    @inbounds for (lhs, i) in s
        for k in 1:nsites(t)
            k_site = readbit(lhs, k)
            rhs = flip(lhs, 1 << (k - 1))
            if haskey(s, rhs)
                dst[i, s[rhs]] = getterm(t, k, k_site)
            end
        end
    end
    return dst
end

function to_matrix!(dst::AbstractMatrix{T}, t::ZTerm, s::Subspace) where T
    @inbounds for (lhs, i) in s
        sigma_z = zero(T)
        for k in 1:nsites(t)
            sigma_z += getterm(t, k, readbit(lhs, k))
        end
        dst[i, i] = sigma_z
    end
    return dst
end

# forward to to_matrix! as fallback
update_term!(H::AbstractMatrix, t::AbstractTerm, s::Subspace) = to_matrix!(H, t, s)

# specialize on sparse matrix
update_term!(H::SparseMatrixCSC, t::AbstractTerm, s::Subspace) = update_term!(H, t, s.subspace_v)

function update_term!(H::SparseMatrixCSC, t::AbstractTerm, subspace_v)
    @inbounds for col in 1:size(H, 1)
        for k in H.colptr[col]:H.colptr[col+1]-1
            row = H.rowval[k]
            lhs = subspace_v[row]
            rhs = subspace_v[col]
            update_nzval!(H.nzval, k, t, col, row, rhs, lhs)
        end
    end
    return H
end

Base.@propagate_inbounds function update_nzval!(nzval, k, t::Hamiltonian, col, row, rhs, lhs)
    for term in t.terms
        update_nzval!(nzval, k, term, col, row, rhs, lhs)
    end
    return nzval
end

Base.@propagate_inbounds function update_nzval!(nzval, k, t::XTerm, col, row, rhs, lhs)
    col == row && return nzval
    mask = lhs ⊻ rhs
    l = log2i(UInt(mask)) + 1
    l_site = lhs & mask
    nzval[k] = getterm(t, l, l_site)
    return nzval
end

Base.@propagate_inbounds function update_nzval!(nzval, k, t::ZTerm, col, row, rhs, lhs)
    col != row && return nzval
    sigma_z = zero(eltype(nzval))
    for k in 1:nsites(t)
        sigma_z += getterm(t, k, readbit(lhs, k))
    end
    nzval[k] = sigma_z
    return nzval
end

Base.@propagate_inbounds getscalarmaybe(x::Vector, k) = x[k]
Base.@propagate_inbounds getscalarmaybe(x::Number, k) = x
Base.@propagate_inbounds getscalarmaybe(x::Nothing, k) = 0
