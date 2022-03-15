function LinearAlgebra.mul!(C::AbstractVecOrMat, A::StepHamiltonian, B::AbstractVecOrMat)
    for (f, term) in zip(A.h.fs, A.h.ts)
        mul!(C, term, B, f(A.t), one(A.t))
    end
    return C
end

# A = sum f_i * A_i
# A * B * α = sum f_i * α * A_i * B
function LinearAlgebra.mul!(C::AbstractVecOrMat, A::StepHamiltonian, B::AbstractVecOrMat, α::Number, β::Number)
    for (f, term) in zip(A.h.fs, A.h.ts)
        mul!(C, term, B, f(A.t) * α, β)
    end
    return C
end
