module RydbergEmulator

using LightGraphs
using LinearAlgebra
using BitBasis
using ExponentialUtilities
using SparseArrays
using Random
using CUDA
import Yao
using Yao: AbstractBlock, AbstractRegister

include("register.jl")
include("measure.jl")
include("unit_disk.jl")
include("hamiltonian.jl")
include("qaoa.jl")
include("qaoa_mis.jl")


@static if CUDA.functional()
    using CUDA.CUSPARSE
    include("cuda.jl")
end

end # module
