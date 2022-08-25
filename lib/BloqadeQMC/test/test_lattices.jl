include("Lattices.jl")
using Plots

# TODO: just playground testing right now
# need more concrete testing

t = 1.
n1 = n2 = 2
PBC = (true, true)

kagome = Kagome(t, n1, n2, PBC)

dij = distance_matrix(kagome; trunc = t, plotting = true)
@show dij