include("Lattices.jl")
using DelimitedFiles

# Define lattice translation vectors
ρ = 1.
a = 4*ρ / sqrt(3)

# a1 must be along x direction
a1 = [a, 0.]
a2 = [a*0.5, a*sqrt(3)*0.5]

# define unit cell basis
r1 = [0., 0.]
r2 = 0.75 * a2
r3 = 0.25 * (a1 + a2)
r4 = 0.5 * a1
r5 = 0.25*a1 + 0.75*a2
r6 = 0.5*a1 + 0.25*a2

r = [r1, r2, r3, r4, r5, r6]

# n1 = number of repitions of the unit cell along a1
# n2 = number of repitions of the unit cell along a2
n1 = 2
n2 = 2

# periodic boundaries along (a1, a2)
PBC = (false, false)

# now define the Custom lattice
L = Custom(a1, a2, n1, n2, r, PBC)
name = "L"
dij = distance_matrix(L)

open("lattice=$(name)_n1=$(n1)_n2=$(n2)_PBC1=$(PBC[1])_PBC2=$(PBC[2]).txt", "w") do io
    writedlm(io, dij)
end

# read in the distance matrix in your QMC as:
# distances = readdlm("lattice=L_n1=2_n2=2_PBC1=false_PBC2=false.txt")
# this gives an Array{Float64,2}, elements are accessed like distances[i,j]