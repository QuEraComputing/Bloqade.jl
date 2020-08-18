abstract type AbstractTerm end
# Ω = 2π * 4 MHz
# C = 
struct RydInteract{T, AtomList <: AbstractVector{<:RydAtom}} <: AbstractTerm
    C::T
    atoms::AtomList
end

const PType = Union{Number, Tuple, Nothing}

struct XTerm{Omega <: PType, Phi <: PType} <: AbstractTerm
    nsites::Int
    Ωs::Omega
    ϕs::Phi
end

struct ZTerm{Delta <: PType} <: AbstractTerm
    nsites::Int
    Δs::Delta
end

struct Hamiltonian{Terms <: Tuple} <: AbstractTerm
    terms::Terms
end

# try to infer number of sites from the input
# this is only necessary for CUDA
to_tuple(xs) = (xs..., ) # make it type stable
to_tuple(xs::Tuple) = xs

XTerm(Ωs::AbstractVector, ϕs::AbstractVector) = XTerm(length(Ωs), to_tuple(Ωs), to_tuple(ϕs))
XTerm(Ωs::Number, ϕs::AbstractVector) = XTerm(length(ϕs), Ωs, to_tuple(ϕs))
XTerm(Ωs::AbstractVector, ϕs::Number) = XTerm(length(Ωs), to_tuple(Ωs), ϕs)

# convenient constructor for simple case
XTerm(n::Int, Ωs::AbstractVector) = XTerm(n, Ωs, nothing)
XTerm(n::Int, Ω::Number) = XTerm(n, Ω, nothing)
XTerm(Ωs::AbstractVector) = XTerm(length(Ωs), to_tuple(Ωs))
ZTerm(Δs::AbstractVector) = ZTerm(length(Δs), to_tuple(Δs))

"""
    nsites(term)

Return the number of sites of given Hamiltonian term.
"""
function nsites end

nsites(t::XTerm) = t.nsites
nsites(t::ZTerm) = t.nsites
nsites(t::Hamiltonian) = nsites(t.terms[1])
nsites(t::RydInteract) = length(t.atoms)

hilbert_space(n::Int) = 0:((1<<n)-1)
hilbert_space(t::AbstractTerm) = hilbert_space(nsites(t))

Base.eltype(t::XTerm) = eltype(t.Ωs)
Base.eltype(t::ZTerm) = eltype(t.Δs)
Base.eltype(t::RydInteract) = typeof(t.C)
Base.eltype(t::Hamiltonian) = promote_type(eltype.(t.terms)...)

Base.getindex(t::Hamiltonian, i::Int) = t.terms[i]

function Base.show(io::IO, t::RydInteract)
    print(io, "C/|r_i - r_j|^6 ")
    printstyled(io, "n_i n_j", color=:light_blue)
end

function Base.show(io::IO, t::XTerm)
    print(io, "Ω ⋅ (e^{iϕ}")
    printstyled(io, "|0)⟨1|", color=:light_blue)
    print(io, "+e^{-iϕ}")
    printstyled(io, "|1⟩⟨0|)", color=:light_blue)
end

function Base.show(io::IO, t::XTerm{<:PType, Nothing})
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

function SparseArrays.SparseMatrixCSC{Tv}(term::AbstractTerm) where {Tv <: Complex}
    N = 1 << nsites(term)
    H = SparseMatrixCOO{Tv}(undef, N, N)
    to_matrix!(H, term)
    return SparseMatrixCSC(H)
end

function SparseArrays.SparseMatrixCSC{Tv}(term::AbstractTerm, s::Subspace) where {Tv <: Complex}
    N = length(s)
    H = SparseMatrixCOO{Tv}(undef, N, N)
    to_matrix!(H, term, s)
    return SparseMatrixCSC(H)
end

SparseArrays.SparseMatrixCSC(term::AbstractTerm, xs...) = SparseMatrixCSC{ComplexF64}(term, xs...)

# full space
# C/|r_i - r_j|^6 n_i n_j
# specialize on COO
function to_matrix!(dst::SparseMatrixCOO, t::RydInteract)
    n = nsites(t)
    @inbounds for i in 1:n, j in 1:i-1
        r_i, r_j = t.atoms[i], t.atoms[j]
        alpha = t.C / distance(r_i, r_j)^6
        # |11⟩⟨11|
        for lhs in itercontrol(nsites(t), [i, j], [1, 1])
            dst[lhs+1, lhs+1] = alpha # this will be accumulated by COO format converter
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

function to_matrix!(dst::SparseMatrixCOO, t::RydInteract, s::Subspace)
    n = nsites(t)
    for (k, lhs) in enumerate(s.subspace_v)
        for i in 1:n, j in 1:i-1
            if (readbit(lhs, i) == 1) && (readbit(lhs, j) == 1)
                r_i, r_j = t.atoms[i], t.atoms[j]
                alpha = t.C / distance(r_i, r_j)^6
                dst[k, k] = alpha
            end
        end
    end
    return dst
end


"""
    update_term!(H, term[, subspace])

Update matrix `H` based on the given Hamiltonian term. This can be faster when the sparse structure of
`H` is known (e.g `H` is a `SparseMatrixCSC`). It fallbacks to `to_matrix!(H, term[, subspace])` if the
sparse structure is unknown.
"""
function update_term! end

# forward to to_matrix! as fallback
update_term!(H::AbstractMatrix, t::AbstractTerm, s::Nothing) = update_term!(H, t)
update_term!(H::SparseMatrixCSC, t::AbstractTerm, s::Nothing) = update_term!(H, t) # disambiguity
update_term!(H::AbstractMatrix, t::AbstractTerm) = to_matrix!(H, t)
update_term!(H::AbstractMatrix, t::AbstractTerm, s::Subspace) = to_matrix!(H, t, s)

# specialize on sparse matrix
update_term!(H::AbstractSparseMatrix, t::AbstractTerm, s::Subspace) = update_term!(H, t, s.subspace_v)

function foreach_nnz(f, H::SparseMatrixCSC)
    for lhs in 1:size(H, 1)
        @inbounds start = H.colptr[lhs]
        @inbounds stop = H.colptr[lhs+1]-1

        for k in start:stop
            @inbounds rhs = H.rowval[k]
            f(k, lhs, rhs)
        end
    end
end

function update_term!(H::AbstractSparseMatrix, t::AbstractTerm)
    nzval = nonzeros(H)
    foreach_nnz(H) do k, col, row
        @inbounds nzval[k] = term_value(t, col-1, row-1, col, row)
    end
    return H
end

function update_term!(H::AbstractSparseMatrix, t::AbstractTerm, subspace_v::AbstractVector)
    nzval = nonzeros(H)
    @inbounds foreach_nnz(H) do k, col, row
        lhs = subspace_v[col]
        rhs = subspace_v[row]

        nzval[k] = term_value(t, lhs, rhs, col, row)
    end
    return H
end

"""
    term_value(term, lhs, rhs, col, row)

Return the value of given term at `H[col, row]` with left basis `lhs` and right basis `rhs`.
For full space, `lhs = col - 1` and `rhs = row - 1`, for subspace, `lhs = subspace_v[col]` and
`rhs = subspace_v[row]`.
"""
function term_value end

@generated function term_value(t::Hamiltonian{Term}, lhs, rhs, col, row) where Term
    ex = Expr(:block)

    push!(ex.args, Expr(:meta, :inline, :propagate_inbounds))
    push!(ex.args, :(val = term_value(t.terms[1], lhs, rhs, col, row)))
    for k in 2:length(Term.parameters)
        push!(ex.args, :(val += term_value(t.terms[$k], lhs, rhs, col, row)))
    end

    push!(ex.args, :val)
    return ex
end

Base.@propagate_inbounds function term_value(t::XTerm, lhs, rhs, col, row)
    col == row && return zero(eltype(t))
    mask = lhs ⊻ rhs
    l = unsafe_log2i(mask) + 1
    l_site = rhs & mask
    return getterm(t, l, l_site)
end

Base.@propagate_inbounds function term_value(t::ZTerm, lhs, rhs, col, row)
    col != row && return zero(eltype(t))
    sigma_z = zero(eltype(t))
    for i in 1:nsites(t)
        sigma_z += getterm(t, i, readbit(lhs, i))
    end
    return sigma_z
end

Base.@propagate_inbounds function term_value(t::RydInteract, lhs, rhs, col, row)
    col != row && return zero(eltype(t))
    # all the nonzeros indices contains
    v = zero(eltype(t))
    n = nsites(t)
    for i in 1:n, j in 1:i-1
        if (readbit(lhs, i) == 1) && (readbit(lhs, j) == 1)
            r_i, r_j = t.atoms[i], t.atoms[j]
            alpha = t.C / distance(r_i, r_j)^6
            v += alpha
        end
    end
    return v
end

Base.@propagate_inbounds getscalarmaybe(x::AbstractVector, k) = x[k]
Base.@propagate_inbounds getscalarmaybe(x::Number, k) = x
Base.@propagate_inbounds getscalarmaybe(x::Tuple, k) = x[k]
Base.@propagate_inbounds getscalarmaybe(x::Nothing, k) = 0

"""
    simple_rydberg(n::Int, ϕ::Number)

Create a simple rydberg hamiltonian that has only [`XTerm`](@ref).
"""
simple_rydberg(n::Int, ϕ::Number) = XTerm(n, one(ϕ), ϕ)

"""
    rydberg_h(C, atoms, Ω, ϕ, Δ)

Create a rydberg hamiltonian

```math
∑ \frac{C}{|r_i - r_j|^6} n_i n_j + Ω σ_x + Δ σ_z
```
"""
function rydberg_h(C, atoms, Ω, ϕ, Δ)
    return RydInteract(C, atoms) + XTerm(length(atoms), Ω, ϕ) + ZTerm(length(atoms), Δ)
end
