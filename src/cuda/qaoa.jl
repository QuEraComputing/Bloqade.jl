function CUDA.cu(Ks::KrylovSubspace)
    KrylovSubspace(Ks.m, Ks.maxiter, Ks.augmented, Ks.beta, cu(Ks.V), cu(Ks.H))
end

function CUDA.cu(qaoa::QAOA)
    cache = CuSparseMatrixCSR(qaoa.cache)
    return QAOA(cu(qaoa.Ks), cache)
end
