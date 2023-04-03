function LinearAlgebra.mul!(C::AbstractVecOrMat, A::StepHamiltonian, B::AbstractVecOrMat)
    fill!(C, zero(eltype(C)))
    for (f, term) in zip(A.h.fs, A.h.ts)
        LinearAlgebra.mul!(C, term, B, f(A.t), one(A.t))
    end
    return C
end

 # default to BloqadeExpression
const backend = @load_preference("backend", "BloqadeExpr")

function LinearAlgebra.mul!(C, A::MultiThreadedMatrix, B, α, β)
    @static if backend == "ParallelMergeCSR"
        ParallelMergeCSR.mul!(C, A.matrix, B, α, β)
    elseif backend == "ThreadedSparseCSR"
        ThreadedSparseCSR.bmul!(C, A.matrix, B, α, β)
    elseif backend == "BloqadeExpr"
        LinearAlgebra.mul!(C, A.matrix, B, α, β)
    else
        throw(ArgumentError("The backend selected is not supported."))
    end
end