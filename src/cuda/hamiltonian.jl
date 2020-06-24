update_term!(H::CuSparseMatrixCSR, t::AbstractTerm, s::Subspace) = update_term!(H, t, s.subspace_v)

function update_term!(dst::CuSparseMatrixCSR, t::AbstractTerm, subspace_v::Vector)
    update_term!(dst, t, CuVector(subspace_v)) # copy to device, if subspace is given on CPU
end

function update_kernel(dst, t::AbstractTerm, subspace_v)
    row = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    if row > length(dst.rowPtr) - 1
        return
    end

    @inbounds for k in dst.rowPtr[row]:dst.rowPtr[row+1]-1
        col = dst.colVal[k]
        lhs = subspace_v[row]
        rhs = subspace_v[col]
        update_nzval!(dst.nzVal, k, t, col, row, rhs, lhs)
    end
end

@generated function update_kernel(dst, t::Hamiltonian{T}, subspace_v) where T
    ex = Expr(:block)
    N = length(T.parameters)
    
    for k in 1:N
        push!(ex.args, :(update_kernel(dst, t.terms[$k], subspace_v)))
    end
    return ex
end

function update_term!(dst::CuSparseMatrixCSR, t::AbstractTerm, subspace_v::CuVector)
    if size(dst, 2) < 256
        threads = size(dst, 2)
        nblocks = 1
    else
        threads = 256
        nblocks = ceil(Int, size(dst, 2) / 256)
    end

    @cuda threads=threads blocks=nblocks update_kernel(dst, t, subspace_v)
    return dst
end
