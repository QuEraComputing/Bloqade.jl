function LinearAlgebra.mul!(C::AbstractVecOrMat, A::StepHamiltonian, B::AbstractVecOrMat)
    fill!(C, zero(eltype(C)))
    for (f, term) in zip(A.h.fs, A.h.ts)
        mul!(C, term, B, f(A.t), one(A.t))
    end
    return C
end
