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
using Base.Threads: nthreads

export Hamiltonian, SumOfLinop
export ThreadedMatrix
export set_backend
export storage_size, to_matrix
export precision_type, highest_type
export add_I, isskewhermitian
#export ValH, get_f # convert StepHamiltonian to ValHamiltonian

include("types.jl")
include("printings.jl")

## things related to low-level linalg of 4 data types:
#  1. PermMatrix
#  2. Diagonal
#  3. SparseMatrixCSC
#  4. SparseMatrixCSR
include("backends/Parallel/linalg_internal.jl")
include("backends/SingleThread/linalg_internal.jl")

include("linalg.jl")
include("preferences.jl")

end 