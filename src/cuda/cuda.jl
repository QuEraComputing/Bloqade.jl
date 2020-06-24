using CUDA
using CUDA.CUSPARSE
using CUDA: CUBLAS
using CUDA: GPUArrays
using CUDA.GPUArrays: AbstractGPUVecOrMat, AbstractGPUArray, AbstractGPUVector

function CUDA.cu(r::RydbergReg{N}) where {N}
    return RydbergReg{N}(cu(r.state), cu(r.subspace))
end

CUDA.cu(s::Subspace) = Subspace(s.map, cu(s.subspace_v))

include("patch.jl")
include("device.jl")
include("hamiltonian.jl")
