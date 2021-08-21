# Copyright 2020 QuEra Computing Inc. All rights reserved.

module CuRydbergEmulator

using RydbergEmulator
using Yao
using CUDA
using CUDA.CUSPARSE
using Adapt
using Reexport
using OrdinaryDiffEq
using SparseArrays
using ContinuousEmulator

using CUDA: CUBLAS
using CUDA: GPUArrays
using CUDA.GPUArrays: AbstractGPUVecOrMat, AbstractGPUArray, AbstractGPUVector
using RydbergEmulator: AbstractTerm
using ContinuousEmulator: EquationCache

@reexport using RydbergEmulator

export cpu
cpu(x) = adapt(Array, x)

include("patch.jl")
include("hamiltonian.jl")
include("ode.jl")

end
