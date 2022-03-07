# write your EaRyd example with Literate.jl here

using EaRyd
using Random

# In this example, we use the Rydberg Emulator to simulate and evolve a fully coherent, 
# strongly interacting system of 9 qubits to observe emergent oscillations in many-body 
# dynamics afer a sudden quench to single-atom resonance. We demonstrate the many-body dynamics 
# with measurements of the domain wall density, which signals the appearance and disappearance of crystalline states.


# We start by building the 1D-Chain 10-atom arrangement, with each atom separated from its neighbor by 5.72 micrometers
# We evaluate the quench dynamics of the Rydberg atom array initially prepared in a product state as the detuning changes to single atom resonance
# After the quence, we observe oscillations of many-body states between the initial and inverted states.

Random.seed!(42)
 # build lattice structure with 10 sites
nsites = 13
atoms = generate_sites(ChainLattice(), nsites, scale=5.72)

# construct Rydberg Hamiltonian with specified Rabi frequency 
h = rydberg_h(atoms;C = 2π * 858386, Ω=4π)

# construct initial product state 
#config = rand(0:1, 9)
#init = product_state(bit"10101010101")
init = zero_state(nsites)

#init = product_state(bit"000000010000000")


# perform discrete time evolution given timestep ts = 0.01 for 120 iterations using Krylov
iteration = 1:120
ts = [0.01 for _ in iteration];
hs = [h for _ in iteration];
prob = KrylovEvolution(init, ts, hs)

# measure observable
clocks = cumsum(ts)
# create empty lists of output expectation values
entropy = zeros(length(iteration)) # entanglement entropy 
domain_mat = zeros(nsites-1, length(iteration)) # domain wall number 
density_mat = zeros(nsites, length(iteration)) # density matrix

corr_mat = zeros(nsites, length(iteration)) # the correlation from the central sites
corr_mat_Z=  zeros(nsites, length(iteration))
kk = Int(floor(nsites/2))+1


for info in prob
    for i in 1:nsites
        density_mat[i, info.step] = expect(put(nsites, i=>Op.n), info.reg)
    end

    for i in 1:nsites-1
        corr = real(expect(put(nsites, (i, i+1)=>kron(Op.n, Op.n)), info.reg))
        obs = density_mat[i, info.step] + density_mat[i+1, info.step] - 2corr
        domain_mat[i, info.step] = obs
    end
    
    for i in 1:nsites
        if i !=kk
        corr= real(expect(put(nsites, (kk, i)=>kron(Op.n, Op.n)), info.reg))
        obs= corr- expect(put(nsites, i=>Op.n), info.reg)* expect(put(nsites, kk=>Op.n), info.reg)
        else
        obs= expect(put(nsites, i=>Op.n), info.reg)- expect(put(nsites, i=>Op.n), info.reg)* expect(put(nsites, kk=>Op.n), info.reg)
        end
        corr_mat[i, info.step]= abs(obs)
    end


    for i in 1:nsites
        if i !=kk
        corr= real(expect(put(nsites, (kk, i)=>kron(Op.Z, Op.Z)), info.reg))
        obs= corr- expect(put(nsites, i=>Op.Z), info.reg)* expect(put(nsites, kk=>Op.Z), info.reg)
        else
        obs= 1- expect(put(nsites, i=>Op.Z), info.reg)* expect(put(nsites, kk=>Op.Z), info.reg)
        end
        corr_mat_Z[i, info.step]= abs(obs)
    end



    rho = density_matrix(info.reg, (1,2,3,4,5))
    entropy[info.step] = von_neumann_entropy(rho)
end

# Plot results 
using CairoMakie
#fig = Figure(size=(10, 5));
#ax = Axis(fig[1, 1])

#for i in 1:nsites
 #   lines!(clocks, density_mat[i, :])
#end
#fig

heatmap(clocks, 1:nsites, (density_mat'))
heatmap(clocks, 1:nsites, corr_mat', c = :solar, size = (1000, 500), colorbar_title = "Rydberg Probability")
heatmap(clocks, 1:nsites, corr_mat_Z', c = :solar, size = (1000, 500), colorbar_title = "Rydberg Probability")



#heatmap(clocks, 1:nsites, domain_mat')
#domain_avg = vec(sum(domain_mat, dims=1)/(nsites-1))
#fig = Figure()
#ax = Axis(fig[1, 1])
#lines!(ax, clocks, domain_avg)
#lines!(ax, clocks, entropy)
#fig
