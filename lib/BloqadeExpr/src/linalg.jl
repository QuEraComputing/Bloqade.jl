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

# Taken directly from ThreadedSparseCSR. ThreadedSparseCSR currently only works with Polyester 0.5, which
# has known issues with thread count > 64. Moving it here allows us to support more recent versions of
# Polyester.
"""
    bmul!(y::AbstractVector, A::SparseMatrixCSR, x::AbstractVector, alpha::Number, beta::Number)
    bmul!(y::AbstractVector, A::SparseMatrixCSR, x::AbstractVector)
Evaluates `y = alpha*A*x + beta*y` (`y = A*x`)
In-place multithreaded version of sparse csr matrix - vector multiplication, using the threading provided by Polyester.jl
"""
function bmul!(y::AbstractVector, A::SparseMatrixCSR, x::AbstractVector, alpha::Number, beta::Number)
    
    A.n == size(x, 1) || throw(DimensionMismatch())
    A.m == size(y, 1) || throw(DimensionMismatch())

    o = getoffset(A)
    
    @batch minbatch = size(y, 1) ÷ matmul_num_threads[] for row in 1:size(y, 1)
        @inbounds begin
            accu = zero(eltype(y))
            for nz in nzrange(A, row)
                col = A.colval[nz] + o
                accu += A.nzval[nz]*x[col]
            end
            y[row] = alpha*accu + beta*y[row]
    
        end
    end

    return y

end

# dispatch to ThreadedSparseCSR
# bmul!(C, A::SparseMatrixCSR, B, α, β) = ThreadedSparseCSR.bmul!(C, A, B, α, β) 
# dispatch to ParallelMergeCSR
bmul!(C, A::Transpose{<:Any, <:SparseMatrixCSC}, B, α, β) = ParallelMergeCSR.mul!(C, A, B, α, β) 