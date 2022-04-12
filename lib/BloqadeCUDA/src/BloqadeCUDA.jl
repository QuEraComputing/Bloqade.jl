# Copyright 2020 QuEra Computing Inc. All rights reserved.

module BloqadeCUDA

using CUDA
using LinearAlgebra
using CUDA.CUSPARSE
using CUDA.CUSPARSE: CuSparseMatrixCSC,
    CuSparseMatrixCSR,
    AbstractCuSparseMatrix

include("patch.jl")

end
