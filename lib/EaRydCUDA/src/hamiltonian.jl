function EaRydKrylovEvolution.storage_size(S::CuSparseMatrixCSR)
    sizeof(S.rowPtr) + sizeof(S.colVal) + sizeof(S.nzVal)
end

function EaRydKrylovEvolution.update_term!(dst::CuSparseMatrixCSR, t::AbstractTerm, s::Subspace)
    update_term!(dst, t, cu(s)) # copy to device, if subspace is given on CPU
end

function thread_layout(H::AbstractCuSparseMatrix)
    return thread_layout(size(H, 2))
end

function thread_layout(len::Int)
    if len < 256
        threads = len
        nblocks = 1
    else
        threads = 256
        nblocks = ceil(Int, len / 256)
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

        H.nzVal[k] = EaRydKrylovEvolution.term_value(t, lhs, rhs, col, row)
    end
    return
end

function update_term_kernel(H, t)
    row = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    @inbounds for k in H.rowPtr[row]:H.rowPtr[row+1]-1
        col = H.colVal[k]

        H.nzVal[k] = EaRydKrylovEvolution.term_value(t, col-1, row-1, col, row)
    end
    return
end

function EaRydKrylovEvolution.update_term!(H::CuSparseMatrixCSR, t::AbstractTerm, ::FullSpace)
    threads, nblocks = thread_layout(H)
    @cuda threads=threads blocks=nblocks update_term_kernel(H, t)
    return H
end

function EaRydKrylovEvolution.update_term!(H::CuSparseMatrixCSR, t::AbstractTerm, s::Subspace{<:CuArray})
    LinearAlgebra.checksquare(H) == length(s) || error("given matrix size does not match subspace size")
    threads, nblocks = thread_layout(H)
    @cuda threads=threads blocks=nblocks update_term_kernel(H, t, s.subspace_v)
    return H
end

function update_dstate_kernel(dstate, state)
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if idx > size(dstate, 1)
        return
    end

    @inbounds begin
        dstate[idx, 1] = state[idx, 2]
        dstate[idx, 2] = -state[idx, 1]
    end
    return
end

function EaRydODEEvolution.update_dstate!(dstate::CuMatrix{<:Real}, state::CuMatrix{<:Real}, ::RealLayout)
    threads, nblocks = thread_layout(size(dstate, 1))
    @cuda threads=threads blocks=nblocks update_dstate_kernel(dstate, state)
    return 
end
