using LuxurySparse
using LinearAlgebra
using SparseArrays


function mul!(Y::Vector, A::PermMatrix, X::Vector, alpha::Number, beta::Number)
    @inbounds for I in axes(X, 1)
        Y[I] = alpha * A.vals[I] * X[A.perm[I]] + beta * Y[I]
    end
    return Y
end

using Yao

mat(put(5, 1=>X))

M = mat(put(20, 1=>Yao.X))

S = SparseMatrixCSC(M)

using BenchmarkTools
st = rand(ComplexF64, 1 << 20)
@benchmark $S * $st
@benchmark $M * $st