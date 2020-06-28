function CUDA.cu(Ks::KrylovSubspace)
end

function CUDA.cu(qaoa::QAOA{N, T}) where {N, T}
    cache = CuSparseMatrixCSR(qaoa.cache)
    Ks = cu(qaoa.Ks)
    expHe = CuVector{T}(undef, Ks.maxiter)
    QAOA{N}(qaoa.ts, cu.(qaoa.term), cache, cu(qaoa.subspace), Ks, expHe)
end
