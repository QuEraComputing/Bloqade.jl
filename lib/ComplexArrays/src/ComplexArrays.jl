module ComplexArrays

using Adapt
using LinearAlgebra

export ComplexArray, ComplexVector, ComplexMatrix

include("types.jl")
include("linalg.jl")
include("broadcast.jl")

end
