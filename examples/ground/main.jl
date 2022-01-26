
# write your EaRyd example with Literate.jl here
# calculate the ground states of a ladder system
using EaRyd
using KrylovKit
using SparseArrays

n_sites= 8
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

#h = rydberg_h(atom_coordinate; C=2.4^6, Δ=4, Ω=1)
C=2.3^6
Δ=2.5 * ones(2*n_sites)
Ω=1* ones(2*n_sites)
Δ[2* n_sites]+= 0.1

h = rydberg_h(atom_coordinate; C=C, Δ=Δ, Ω=Ω)

H0 = SparseMatrixCSC(h)
vals, vecs, info = KrylovKit.eigsolve(H0,  4, :SR)
state1= ArrayReg(vecs[1]) 

output_mat1= zeros(Float64, 2*n_sites)

for i in 1:2*n_sites
    
    output_mat1[i] = real(expect(put(2*n_sites, i=>Op.n), state1))

end

println(output_mat1)

n_diff= output_mat1[1:n_sites]- output_mat1[n_sites+1: 2*n_sites]
plot(1: n_sites, n_diff)

# next step, plot things like Jin 
