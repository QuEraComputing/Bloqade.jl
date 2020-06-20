module RydbergEmulator

using Printf
using BitBasis
using SparseArrays

export RydInteract, RydAtom, XTerm, ZTerm, Hamiltonian

include("atoms.jl")
include("hamiltonian2.jl")

# using LightGraphs
# using LinearAlgebra
# using BitBasis
# using ExponentialUtilities
# using SparseArrays
# using Random
# using CUDA
# import Yao
# using Yao: AbstractBlock, AbstractRegister

# include("register.jl")
# include("measure.jl")
# include("unit_disk.jl")
# include("hamiltonian.jl")
# include("qaoa.jl")
# include("qaoa_mis.jl")


# @static if CUDA.functional()
#     using CUDA.CUSPARSE
#     include("cuda.jl")
# end

end # module
