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
using LinearAlgebra

using CUDA: CUBLAS
using CUDA: GPUArrays
using CUDA.GPUArrays: AbstractGPUVecOrMat, AbstractGPUArray, AbstractGPUVector
using RydbergEmulator: AbstractTerm, KrylovEmulationCache, PrecisionAdaptor
using ContinuousEmulator: EquationCache, ContinuousOptions

@reexport using RydbergEmulator
@reexport using ContinuousEmulator

export cpu
cpu(x) = adapt(Array, x)

LinearAlgebra.normalize!(r::RydbergReg{N, 1, <:CuArray}) where N = normalize!(vec(r.state))

Adapt.adapt_storage(::PrecisionAdaptor{P}, x::CuArray{<:Real}) where P = convert(CuArray{P}, x)
Adapt.adapt_storage(::PrecisionAdaptor{P}, x::CuArray{<:Complex}) where P = convert(CuArray{Complex{P}}, x)

include("hamiltonian.jl")
include("ode.jl")

end
