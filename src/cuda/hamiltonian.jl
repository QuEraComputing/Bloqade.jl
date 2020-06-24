function update_term!(dst::CuSparseMatrixCSR, n::Int, subspace_v::CuVector)
end

function update_hamiltonian!(dst::CuSparseMatrixCSR, n::Int, subspace_v::CuVector, Ω, ϕ)
    function kernel(rowPtr, colVal, nzVal, subspace_v, Ω, ϕ)
        row = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        if row > length(rowPtr) - 1
            return
        end

        for k in rowPtr[row]:rowPtr[row+1]-1
            col = colVal[k]
            lhs = subspace_v[row]
            rhs = subspace_v[col]

            mask = lhs ⊻ rhs
            if lhs & mask == 0
                nzVal[k] = Ω * exp(im * ϕ)
            else
                nzVal[k] = Ω * exp(-im * ϕ)
            end
            # update_x_term!(nzVal, lhs, rhs, k, Ω, ϕ)
        end
        return
    end
    
    if size(dst, 2) < 256
        threads = size(dst, 2)
        nblocks = 1
    else
        threads = 256
        nblocks = ceil(Int, size(dst, 2) / 256)
    end

    @cuda threads=threads blocks=nblocks kernel(dst.rowPtr, dst.colVal, dst.nzVal, subspace_v, Ω, ϕ)
    return dst
end

function update_hamiltonian!(dst::CuSparseMatrixCSR, n::Int, subspace_v::CuVector, Ω, ϕ, Δ)
    function kernel(rowPtr, colVal, nzVal, subspace_v, Ω, ϕ, Δ)
        row = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        if row > length(rowPtr) - 1
            return
        end

        for k in rowPtr[row]:rowPtr[row+1]-1
            col = colVal[k]
            lhs = subspace_v[row]
            rhs = subspace_v[col]

            if row == col
                update_z_term!(nzVal, lhs, n, k, Δ)
            else
                update_x_term!(nzVal, lhs, rhs, k, Ω, ϕ)
            end
        end
        return
    end

    if size(dst, 2) < 256
        threads = size(dst, 2)
        nblocks = 1
    else
        threads = 256
        nblocks = ceil(Int, size(dst, 2) / 256)
    end

    @cuda threads=threads blocks=nblocks kernel(dst.rowPtr, dst.colVal, dst.nzVal, subspace_v, Ω, ϕ, Δ)
    return dst
end
