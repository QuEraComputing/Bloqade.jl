
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

h = rydberg_h(atom_coordinate; C=2.4^6, Δ=3.5, Ω=1)

H0 = SparseMatrixCSC(h)
vals, vecs, info = KrylovKit.eigsolve(H0,  4, :SR)
state1= ArrayReg(vecs[1]) 
state2= ArrayReg(vecs[2])

state3= ArrayReg((vecs[1]-vecs[2])/sqrt(2))


output_mat1= zeros(Float64, 2*n_sites)
output_mat2= zeros(Float64, 2*n_sites)
output_mat3= zeros(Float64, 2*n_sites)



for i in 1:2*n_sites
    observe1 = real(expect(put(2*n_sites, i=>Op.n), state1))
    observe2 = real(expect(put(2*n_sites, i=>Op.n), state2))
    observe3 = real(expect(put(2*n_sites, i=>Op.n), state3))

    output_mat1[i] = observe1
    output_mat2[i] = observe2
    output_mat3[i] = observe3

end

println(output_mat1)
println(output_mat2)
println(output_mat3)
