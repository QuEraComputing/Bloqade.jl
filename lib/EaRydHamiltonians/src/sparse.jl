"""
    getterm(terms, k, k_site)

Get the value of k-th local term in `terms`
given the site configuration as `k_site`.
"""
function getterm end

getterm(t::Negative, k, k_site) = -getterm(t.term, k, k_site)

function getterm(t::XTerm, k, k_site)
    if k_site == 0
        return getscalarmaybe(t.Ωs, k)/2 * exp(im * getscalarmaybe(t.ϕs, k))
    else
        return getscalarmaybe(t.Ωs, k)/2 * exp(-im * getscalarmaybe(t.ϕs, k))
    end
end

function getterm(t::XTerm{<:Any, Nothing}, k, k_site)
    return getscalarmaybe(t.Ωs, k)/2
end

function getterm(t::ZTerm, k, k_site)
    if k_site == 0
        return getscalarmaybe(t.Δs, k)
    else
        return -getscalarmaybe(t.Δs, k)
    end
end

function getterm(t::NTerm, k, k_site)
    if k_site == 1
        return getscalarmaybe(t.Δs, k)
    else
        return 0
    end
end

function getterm(t::Hamiltonian, k, k_site)
    error("composite Hamiltonian term cannot be indexed")
end

space_size(term::AbstractTerm, s::FullSpace) = 1 << nsites(term)
space_size(::AbstractTerm, s::Subspace) = length(s)

function SparseArrays.SparseMatrixCSC{Tv, Ti}(term::AbstractTerm, s::AbstractSpace=fullspace) where {Tv, Ti}
    N = space_size(term, s)
    colptr, rowval = sparse_skeleton_csc(Ti, term, s)
    H = SparseMatrixCSC{Tv, Ti}(N, N, colptr, rowval, Vector{Tv}(undef, length(rowval)))
    update_term!(H, term, s)
    return H
end

function SparseArrays.SparseMatrixCSC{Tv, Ti}(term::Hamiltonian, s::AbstractSpace=fullspace) where {Tv, Ti}
    return sum(SparseMatrixCSC{Tv, Ti}(t, s) for t in term.terms)
end

function SparseArrays.SparseMatrixCSC{Tv, Ti}(t::Negative, s::AbstractSpace=fullspace) where {Tv, Ti}
    H = SparseMatrixCSC{Tv, Ti}(t.term, s)
    lmul!(-one(Tv), H.nzval)
    return H
end

# NOTE: we use BlasInt as Ti, since MKLSparse uses BlasInt as Ti
# on some devices BlasInt != Int, thus this is necessary to trigger
# dispatch to MKLSparse
SparseArrays.SparseMatrixCSC{Tv}(term::AbstractTerm, s::AbstractSpace=fullspace) where Tv = SparseMatrixCSC{Tv, BlasInt}(term, s)
SparseArrays.SparseMatrixCSC(term::AbstractTerm, s::AbstractSpace=fullspace) = SparseMatrixCSC{ComplexF64}(term, s)

sparse_skeleton_csc(t::AbstractTerm, s::AbstractSpace=fullspace) = sparse_skeleton_csc(Int, t, s)

function sparse_skeleton_csc(::Type{Ti}, t::XTerm, s::FullSpace) where Ti
    natoms = nsites(t)
    sz = hilbert_space(natoms)
    colptr = Ti[natoms*i+1 for i in 0:length(sz)]
    rowval = Vector{Ti}(undef, natoms * length(sz))
    @inbounds @batch for col in sz
        colrowval = map(1:natoms) do k
            rhs = flip(col, 1 << (k - 1))
            return rhs + 1
        end
        rowval[natoms*col+1:natoms*col+natoms] = sort!(colrowval)
    end
    return colptr, rowval
end

function sparse_skeleton_csc(::Type{Ti}, t::XTerm, s::Subspace) where Ti
    colptr = Vector{Ti}(undef, length(s)+1)
    rowval_list = Vector{Vector{Ti}}(undef, length(s))
    colptr[1] = 1

    @inbounds @batch for i in 1:length(s)
        lhs = s.subspace_v[i]
        colrowval = Ti[]
        for k in 1:nsites(t)
            rhs = flip(lhs, 1 << (k - 1))
            if haskey(s, rhs)
                push!(colrowval, s[rhs])
            end
        end
        sort!(colrowval)
        rowval_list[i] = colrowval
    end

    rowval = Vector{Ti}(undef, sum(length, rowval_list))
    ptr = 1
    for col in 1:length(rowval_list) # O(length(s))
        rowval_col = rowval_list[col]
        rowval_nnz = length(rowval_col)
        rowval[ptr:ptr+rowval_nnz-1] = rowval_col
        colptr[col+1] = colptr[col] + rowval_nnz
        ptr += rowval_nnz
    end
    return colptr, rowval
end

# Diagonal
function sparse_skeleton_csc(::Type{Ti}, t::Union{RydInteract, ZTerm, NTerm}, s::AbstractSpace) where Ti
    m = space_size(t, s)
    colptr = Vector{Ti}(1:(m+1))
    rowval = Vector{Ti}(1:m)
    return colptr, rowval
end

"""
    update_term!(H, term[, space=fullspace])

Update matrix `H` based on the given Hamiltonian term. This can be faster when the sparse structure of
`H` is known (e.g `H` is a `SparseMatrixCSC`).
"""
function update_term! end

Base.@propagate_inbounds function foreach_nnz(f, H::SparseMatrixCSC)
    @batch for lhs in 1:size(H, 1)
        @inbounds start = H.colptr[lhs]
        @inbounds stop = H.colptr[lhs+1]-1

        for k in start:stop
            @inbounds rhs = H.rowval[k]
            f(k, lhs, rhs)
        end
    end
end

function update_term!(H::AbstractSparseMatrix, t::AbstractTerm, ::FullSpace=fullspace)
    nzval = nonzeros(H)
    foreach_nnz(H) do k, col, row
        @inbounds nzval[k] = term_value(t, col-1, row-1, col, row)
    end
    return H
end

function update_term!(H::AbstractSparseMatrix, t::AbstractTerm, s::Subspace)
    nzval = nonzeros(H)
    @inbounds foreach_nnz(H) do k, col, row
        lhs = s.subspace_v[col]
        rhs = s.subspace_v[row]
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

Base.@propagate_inbounds function term_value(t::OrNegative{<:XTerm}, lhs, rhs, col, row)
    col == row && return zero(eltype(t))
    mask = lhs ⊻ rhs
    l = unsafe_log2i(mask) + 1
    l_site = rhs & mask
    return getterm(t, l, l_site)
end

Base.@propagate_inbounds function term_value(t::OrNegative{<:Union{ZTerm, NTerm}}, lhs, rhs, col, row)
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
Base.@propagate_inbounds getscalarmaybe(x::Tuple, k) = x[k]
@inline getscalarmaybe(x::Number, k) = x
@inline getscalarmaybe(x::Nothing, k) = 0
@inline getscalarmaybe(x, k) = x

"""
    simple_rydberg(n::Int, ϕ::Number)

Create a simple rydberg hamiltonian that has only [`XTerm`](@ref).
"""
simple_rydberg(n::Int, ϕ::Number) = XTerm(n, one(ϕ), ϕ)

"""
    rydberg_h(atoms, [C=2π * 109.133 * MHz*µm^6], Ω, ϕ, Δ)

Create a rydberg hamiltonian, shorthand for
`RydInteract(C, atoms) + XTerm(length(atoms), Ω, ϕ) + ZTerm(length(atoms), Δ)`

```math
∑ \\frac{C}{|r_i - r_j|^6} n_i n_j + Ω σ_x - Δ σ_n
```
"""
function rydberg_h(atoms, C, Ω, ϕ, Δ)
    return RydInteract(atoms, C) + XTerm(length(atoms), Ω, ϕ) - NTerm(length(atoms), Δ)
end

function rydberg_h(atoms, Ω, ϕ, Δ)
    return RydInteract(atoms) + XTerm(length(atoms), Ω, ϕ) - NTerm(length(atoms), Δ)
end

function rydberg_h(atoms; Ω, Δ, C=nothing, ϕ=nothing)
    term = - NTerm(length(atoms), Δ)
    if C === nothing
        term += RydInteract(atoms)
    else
        term += RydInteract(atoms, C)
    end

    if ϕ === nothing
        term += XTerm(length(atoms), Ω)
    else
        term += XTerm(length(atoms), Ω, ϕ)
    end
    return term
end

function is_time_dependent(h::XTerm)
    return !(h.Ωs isa ConstParamType || isnothing(h.Ωs)) ||
        !(h.ϕs isa ConstParamType || isnothing(h.ϕs))
end

function is_time_dependent(h::Union{ZTerm, NTerm})
    return !(h.Δs isa ConstParamType)
end

# NOTE: we currently assume atom positions are constant
# might change in the future
is_time_dependent(h::RydInteract) = false
is_time_dependent(t::Negative) = is_time_dependent(t.term)
function is_time_dependent(h::Hamiltonian)
    any(is_time_dependent, h.terms)
end

# Time dependent Term
attime(x::Nothing, t::Real) = nothing
attime(x::Number, t::Real) = x
attime(x, t::Real) = x(t)
attime(x::AbstractArray, t::Real) = attime.(x, t)
attime(x::NTuple{N, <:Number}, t::Real) where N = attime.(x, t)

function (tm::XTerm)(t::Real)
    return XTerm(tm.nsites, attime(tm.Ωs, t), attime(tm.ϕs, t))
end

function (tm::ZTerm)(t::Real)
    return ZTerm(tm.nsites, attime(tm.Δs, t))
end

function (tm::NTerm)(t::Real)
    return NTerm(tm.nsites, attime(tm.Δs, t))
end

(tm::Negative)(t::Real) = Negative(tm.term(t))

function (tm::Hamiltonian)(t::Real)
    return Hamiltonian(map(x->x(t), tm.terms))
end

# fallback to constants
(tm::AbstractTerm)(t::Real) = tm
