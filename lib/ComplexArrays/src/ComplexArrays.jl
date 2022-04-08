module ComplexArrays

using Adapt
using LinearAlgebra

export ComplexArray, ComplexVector, ComplexMatrix

const ComplexVector{T} = ComplexArray{T, 1}
const ComplexMatrix{T} = ComplexArray{T, 2}

include("types.jl")

end
