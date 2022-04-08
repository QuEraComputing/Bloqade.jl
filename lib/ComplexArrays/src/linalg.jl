function LinearAlgebra.mul!(y::ComplexVector, A::AbstractVecOrMat, x::ComplexVector, alpha::Number, beta::Number)
    mul!(y.storage, A, x.storage, alpha, beta)
    return y
end
