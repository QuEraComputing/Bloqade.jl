function update_term!(dst::CuSparseMatrixCSR, t::AbstractTerm, subspace_v::Vector)
    update_term!(dst, t, CuVector(subspace_v)) # copy to device, if subspace is given on CPU
end

function update_term!(dst::CuSparseMatrixCSR, t::AbstractTerm, subspace_v::CuVector)
    function kernel(rowPtr, colVal, nzVal, subspace_v)
        row = (blockIdx().x - 1) * blockDim().x + threadIdx().x

        if row > length(rowPtr) - 1
            return
        end

        @inbounds for k in rowPtr[row]:rowPtr[row+1]-1
            col = colVal[k]
            lhs = subspace_v[row]
            rhs = subspace_v[col]
            update_nzval!(nzVal, k, t, col, row, rhs, lhs)
        end
    end

    if size(dst, 2) < 256
        threads = size(dst, 2)
        nblocks = 1
    else
        threads = 256
        nblocks = ceil(Int, size(dst, 2) / 256)
    end

    @cuda threads=threads blocks=nblocks kernel(dst.rowPtr, dst.colVal, dst.nzVal, subspace_v)
    return dst
end
