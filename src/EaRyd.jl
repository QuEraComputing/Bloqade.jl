module EaRyd

export FangliLattice, LatticeRow, site_positions

using Reexport

@reexport using RydbergEmulator
# @reexport using ContinuousEmulator
# @reexport using Measurements: ±, Measurement

# using CUDA
# @static if CUDA.functional()
#     @reexport using CuRydbergEmulator
# end

include("fangli.jl")

end
