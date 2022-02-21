function SparseArrays.SparseMatrixCSC{Tv, Ti}(term::AbstractTerm, s::AbstractSpace=fullspace) where {Tv, Ti}
    N = space_size(term, s)
    colptr, rowval = sparse_skeleton_csc(Ti, term, s)
    H = SparseMatrixCSC{Tv, Ti}(N, N, colptr, rowval, Vector{Tv}(undef, length(rowval)))
    update_term!(H, term, s)
    return H
end

# function SparseArrays.SparseMatrixCSC{Tv, Ti}(term::Hamiltonian, s::AbstractSpace=fullspace) where {Tv, Ti}
#     return sum(SparseMatrixCSC{Tv, Ti}(t, s) for t in term.terms)
# end

function SparseArrays.SparseMatrixCSC{Tv, Ti}(t::Negative, s::AbstractSpace=fullspace) where {Tv, Ti}
    H = SparseMatrixCSC{Tv, Ti}(t.term, s)
    lmul!(-one(Tv), H.nzval)
    return H
end

# NOTE: we use BlasInt as Ti, since MKLSparse uses BlasInt as Ti
# on some devices BlasInt != Int, thus this is necessary to trigger
# dispatch to MKLSparse
SparseArrays.SparseMatrixCSC{Tv}(term::AbstractTerm, s::AbstractSpace=fullspace) where Tv = SparseMatrixCSC{Tv, BlasInt}(term, s)
SparseArrays.SparseMatrixCSC(term::AbstractTerm, s::AbstractSpace=fullspace) = SparseMatrixCSC{Complex{eltype(term)}}(term, s)

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