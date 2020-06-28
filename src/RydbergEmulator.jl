module RydbergEmulator

using Printf
using BitBasis
using SparseArrays
using LuxurySparse
using LightGraphs
using LinearAlgebra
using OrderedCollections
using ExponentialUtilities
using Random
using CUDA
import Yao

using LinearAlgebra: BlasReal, BlasComplex

export RydInteract, RydAtom, XTerm, ZTerm, Hamiltonian, QAOA
export to_matrix, to_matrix!, update_term!, simple_rydberg, rydberg_h

include("atoms.jl")
include("subspace.jl")
include("hamiltonian.jl")

include("register.jl")
include("measure.jl")

include("unit_disk.jl")

include("qaoa2.jl")

# using LightGraphs
# using LinearAlgebra
# using BitBasis
# using ExponentialUtilities
# using SparseArrays
# using Random
# using CUDA
# import Yao
# using Yao: AbstractBlock, AbstractRegister

# include("qaoa.jl")
# include("qaoa_mis.jl")

@static if CUDA.functional()
    include("cuda/cuda.jl")
end

end # module
