module Lowlevel

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

## things related to low-level linalg of three types:
#  1. PermMatrix
#  2. Diagonal
#  3. SparseMatrixCSC
#  4. SparseMatrixCSR
include("backends/Parallel/linalg_internal.jl")
include("backends/SingleThread/linalg_internal.jl")

include("linalg.jl")
include("preferences.jl")

end 