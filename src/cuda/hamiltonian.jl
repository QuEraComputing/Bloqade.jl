function update_term!(dst::CuSparseMatrixCSR, t::AbstractTerm, subspace_v::Vector)
    update_term!(dst, t, CuVector(subspace_v)) # copy to device, if subspace is given on CPU
end

function thread_layout(H)
    if size(H, 2) < 256
        threads = size(H, 2)
        nblocks = 1
    else
        threads = 256
        nblocks = ceil(Int, size(H, 2) / 256)
    end
    return threads, nblocks
end

function update_term_kernel(H, t, subspace_v)
    row = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    if row > length(H.rowPtr) - 1
        return
    end

    @inbounds for k in H.rowPtr[row]:H.rowPtr[row+1]-1
        col = H.colVal[k]
        lhs = subspace_v[col]
        rhs = subspace_v[row]

        H.nzVal[k] = term_value(t, lhs, rhs, col, row)
    end
    return
end

function update_term_kernel(H, t)
    row = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    if row > length(H.rowPtr) - 1
        return
    end

    @inbounds for k in H.rowPtr[row]:H.rowPtr[row+1]-1
        col = H.colVal[k]
        @inbounds H.nzVal[k] = term_value(t, col-1, row-1, col, row)
    end
    return
end

function update_term!(H::CuSparseMatrixCSR, t::AbstractTerm)
    threads, nblocks = thread_layout(H)
    @cuda threads=threads blocks=nblocks update_term_kernel(H, t)
    return H
end

function update_term!(H::CuSparseMatrixCSR, t::AbstractTerm, subspace_v::CuVector)
    threads, nblocks = thread_layout(H)
    @cuda threads=threads blocks=nblocks update_term_kernel(H, t, subspace_v)
    return H
end
