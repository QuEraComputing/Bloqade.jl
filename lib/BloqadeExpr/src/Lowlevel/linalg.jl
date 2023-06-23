#**************************************************************
#  Here, one can find the linalg API for  
#  1. Hamiltonian/SumOfLinop
#  2. Binding of standard LinearAlgebra functions 
#     to the backend of choice for ThreadedMatrix
#**************************************************************


## mul!
#--------------------------------
function LinearAlgebra.mul!(C::AbstractVecOrMat, A::SumOfLinop, B::AbstractVecOrMat)
    fill!(C, zero(eltype(C)))
    for (f, term) in zip(A.fvals, A.h.ts)
        mul!(C, term, B, f, one(f))
    end
    return C
end

## additionals, maybe we don't need this.
function Base.:*(a::Number, b::SumOfLinop)
    return SumOfLinop(b.fvals .* a, b.h)
end
Base.:*(n, m::T) where {T <: ThreadedMatrix} = n * m.matrix

function Base.:+(a::SumOfLinop, b::SumOfLinop)
    if !(a === b)
        error("two SumOfLinop must share the same static terms ")
    end
    return SumOfLinop(a.fvals + b.fvals, a.h)
end

function Base.:-(a::SumOfLinop, b::SumOfLinop)
    if !(a === b)
        error("two SumOfLinop must share the same static terms ")
    end
    return SumOfLinop(a.fvals - b.fvals, a.h)
end



# if BloqadeExpr was the backend of choice, then ThreadedMatrix will have the SparseMatrixCSC type
# in which case we just dispatch to standard mul!
LinearAlgebra.mul!(C, A::ThreadedMatrix{<:SparseMatrixCSC}, B, α, β) = SparseArrays.mul!(C, A.matrix, B, α, β) 

# if PMCSR/TSCSR was the backend of choice, then direct to bmul!'s below
LinearAlgebra.mul!(C, A::ThreadedMatrix, B, α, β) = bmul!(C, A.matrix, B, α, β)

##-------------------------------- mul!

## opnorm()
# --------------------------------
function LinearAlgebra.opnorm(h::SumOfLinop, p = 2)
    return opnorm(to_matrix(h), p)
end

##---------------------------------


## tr()
# --------------------------------
function LinearAlgebra.tr(A::SumOfLinop)
    return sum(zip(A.fvals, A.h.ts)) do (f, t)
        return f * tr(t)
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