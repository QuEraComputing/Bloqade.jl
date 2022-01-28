# write your EaRyd example with Literate.jl here
# calculate the ground states of a ladder system
using EaRyd
using KrylovKit
using SparseArrays
using Plots

# we create a two ladder of atoms with 8 sites in each row
n_sites= 8

# x coordinate for each row
ax=1
L = (n_sites - 1)* ax
atomx= collect(0.0:ax:L)

# y coordinate for the two rows
ay1=0
ay2= 2

# generate the atom sites one by one, since we dont have rectangular lattice now
atom_coordinate=[(0.0, 0.0)]
for ii= 2: n_sites
    push!(atom_coordinate, (atomx[ii], ay1))
end
for ii= 1: n_sites
    push!(atom_coordinate, (atomx[ii], ay2))
end

# Hamiltonian parameters
C=2.3^6      
Δ=4 * ones(2*n_sites)     
Ω=1* ones(2*n_sites)
Δ[2* n_sites]+= -0.05  # addtional detuning term to break the symmetry 

# create the Hamiltonian 
h = rydberg_h(atom_coordinate; C=C, Δ=Δ, Ω=Ω)

# the sparse matrix for the hamiltonian 
H0 = SparseMatrixCSC(h)

# diagonize the hamiltonian 
vals, vecs, info = KrylovKit.eigsolve(H0,  4, :SR)
state1= ArrayReg(vecs[1]) 

# measruing the magnetization for the ground state 
output_mat1= zeros(Float64, 2*n_sites)
for i in 1:2*n_sites
     output_mat1[i] = real(expect(put(2*n_sites, i=>Op.n), state1))
end
println(output_mat1)

# plot the difference of the magnetization in the two rows 
n_diff= output_mat1[1:n_sites]- output_mat1[n_sites+1: 2*n_sites]
plot(1: n_sites, n_diff)

