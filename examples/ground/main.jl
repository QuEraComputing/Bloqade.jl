
# write your EaRyd example with Literate.jl here


# calculate the ground states of a ladder system
using EaRyd
using KrylovKit
using SparseArrays

n_sites= 10
L = (n_sites - 1)* ax
atomx= collect(0.0:ax:L)


ay1=0
ay2= 2

atom_coordinate=[(0.0, 0.0)]

for ii= 2: n_sites
    push!(atom_coordinate, (atomx[ii], ay1))
end

for ii= 1: n_sites
    push!(atom_coordinate, (atomx[ii], ay2))
end

h = rydberg_h(atom_coordinate; C=1.5^6, Δ=3, Ω=1)

H0 = SparseMatrixCSC(h)
vals, vecs, info = KrylovKit.eigsolve(H0,  1, :SR)
state1= ArrayReg(vecs[1]) 

output_mat1= zeros(Float64, 2*n_sites)


for i in 1:2*n_sites
    observe = real(expect(put(2*n_sites, i=>Op.n), state1))
    output_mat1[i] = observe
end

print(output_mat1)