module RydbergEmulator

using Printf
using BitBasis
using SparseArrays
using LightGraphs
using LinearAlgebra
using OrderedCollections
using CUDA

using LinearAlgebra: BlasReal, BlasComplex

export RydInteract, RydAtom, XTerm, ZTerm, Hamiltonian
export to_matrix, to_matrix!, update_term!

include("atoms.jl")
include("subspace.jl")
include("hamiltonian.jl")

include("register.jl")
include("measure.jl")

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


@static if CUDA.functional()
    include("cuda/cuda.jl")
end

end # module
