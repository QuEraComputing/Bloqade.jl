# Copyright 2020 QuEra Computing Inc. All rights reserved.

module CuRydbergEmulator

using RydbergEmulator
using CUDA
using CUDA.CUSPARSE
using CUDA: CUBLAS
using CUDA: GPUArrays
using CUDA.GPUArrays: AbstractGPUVecOrMat, AbstractGPUArray, AbstractGPUVector
using Adapt
using Reexport
using RydbergEmulator: AbstractTerm

@reexport using RydbergEmulator

export cpu
cpu(x) = adapt(Array, x)

include("patch.jl")
include("hamiltonian.jl")

end
