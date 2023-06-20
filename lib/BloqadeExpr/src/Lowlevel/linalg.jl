#**************************************************************
#  Here, one can find the linalg API for  
#  1. Hamiltonian/StepHamiltonian
#  2. Binding of standard LinearAlgebra functions to the backend of choice
#**************************************************************


## mul!
#--------------------------------
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

##-------------------------------- mul!


## tr()
# --------------------------------
function LinearAlgebra.tr(A::StepHamiltonian)
    return sum(zip(A.h.fs, A.h.ts)) do (f, t)
        return f(A.t) * tr(t)
    end
end

# [TODO] parallel trace (btrace)
# if BloqadeExpr was the backend of choice, then ThreadedMatrix will have the SparseMatrixCSC type
# in which case we just dispatch to single threaded tr
##LinearAlgebra.tr(C, A::ThreadedMatrix{<:SparseMatrixCSC}, B, α, β) = trace(A.matrix) 

# if PMCSR/TSCSR was the backend of choice, then direct to bmul!'s below
##LinearAlgebra.tr(C, A::ThreadedMatrix, B, α, β) = btrace(A.matrix)

## [NOTE] currently only single thread for trace is implemented, so we use the exposed LinearAlgebra.tr() 
#  that bind to trace()
LinearAlgebra.tr(A::ThreadedMatrix) = tr(A.matrix)

##--------------------------------  tr()