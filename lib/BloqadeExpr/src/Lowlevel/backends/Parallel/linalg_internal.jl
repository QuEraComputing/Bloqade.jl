
# Parallelized Permutation Matrix multiplication
function bmul!(C::AbstractVector, A::PermMatrix, B::AbstractVector, α::Number, β::Number)

    A_vals = A.vals
    A_perm = A.perm
    @batch for i in eachindex(B)
        perm_mul!(C, A_vals, A_perm, B, α, β, i)
    end
    return C
end
@inline function perm_mul!(C, A_vals, A_perm, B, α, β, i)
    @inbounds C[i] = A_vals[i] * B[A_perm[i]] * α + β * C[i]
end

# Parallelized Diagonal multiplication
function bmul!(C, A::Diagonal, B, α, β)
    d = A.diag
    @batch for i in eachindex(B)
        @inbounds C[i] = d[i] * B[i] * α + β * C[i] 
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
    
    @batch minbatch = size(y, 1) ÷ nthreads() for row in 1:size(y, 1)
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
#bmul!(C, A::Transpose{<:Any, <:SparseMatrixCSC}, B, α, β) = ParallelMergeCSR.mul!(C, A, B, α, β) 
function bmul!(C, A::Transpose{<:Any, <:SparseMatrixCSC}, B, α, β)
    println("calling PMCSR.mul!")
    return ParallelMergeCSR.mul!(C, A, B, α, β) 
end