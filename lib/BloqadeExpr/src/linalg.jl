function LinearAlgebra.mul!(C::AbstractVecOrMat, A::StepHamiltonian, B::AbstractVecOrMat)
    fill!(C, zero(eltype(C)))
    for (f, term) in zip(A.h.fs, A.h.ts)
        mul!(C, term, B, f(A.t), one(A.t))
    end
    return C
end

# if BloqadeExpr was the backend of choice, then ThreadedMatrix will have the SparseMatrixCSC type
# in which case we just dispatch to standard mul!
LinearAlgebra.mul!(C, A::ThreadedMatrix{<:SparseMatrixCSC}, B, α, β) = SparseArrays.mul!(C, A.matrix, B, α, β) 
# if PMCSR/TSCSR was the backend of choice, then direct to bmul!'s below
LinearAlgebra.mul!(C, A::ThreadedMatrix, B, α, β) = bmul!(C, A.matrix, B, α, β)

# Parallelized Permutation Matrix multiplication
function bmul!(Y::AbstractVector, A::PermMatrix, X::AbstractVector, alpha::Number, beta::Number)

    length(X) == size(A, 2) || throw(DimensionMismatch("input X length does not match PermMatrix A"))
    length(Y) == size(A, 2) || throw(DimensionMismatch("output Y length does not match PermMatrix A"))

    @inbounds @batch for I in eachindex(X)
        Y[I] = A.vals[I] * X[A.perm[I]] * alpha + beta * Y[I]
    end
    return Y
end

# Parallelized Diagonal multiplication
function bmul!(C, A::Diagonal, B, α, β)
    @batch for i in eachindex(B)
        @inbounds C[i] = A.diag[i] * B[i] * α + β * C[i] 
    end
    return C
end

# dispatch to ThreadedSparseCSR
bmul!(C, A::SparseMatrixCSR, B, α, β) = ThreadedSparseCSR.bmul!(C, A, B, α, β) 
# dispatch to ParallelMergeCSR
bmul!(C, A::Transpose{<:Any, <:SparseMatrixCSC}, B, α, β) = ParallelMergeCSR.mul!(C, A, B, α, β) 