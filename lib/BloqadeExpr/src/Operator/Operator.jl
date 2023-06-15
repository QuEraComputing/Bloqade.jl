module Operator

    using SparseArrays
    using LuxurySparse
    using SparseMatricesCSR
    using ParallelMergeCSR
    using Polyester
    using Preferences
    using Adapt
    using LaTeXStrings
    using LinearAlgebra

    export Hamiltonian, StepHamiltonian
    export ThreadedMatrix
    export set_backend
    export storage_size, to_matrix

    include("types.jl")
    include("printings.jl")
    include("linalg.jl")
    include("preferences.jl")

end 