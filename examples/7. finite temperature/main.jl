# # Finite Temperature ED
# ## Background

# In this short example, we show how to use exact diagonalization to compute properties of thermal states for small system sizes. For ground state properties, please see the tutorial on Adiabatic Evolution.

# Let's start by importing the required libraries:

using Bloqade: rydberg_density
using BloqadeLattices: generate_sites, ChainLattice
using BloqadeExpr: rydberg_h 

using KrylovKit
using LinearAlgebra
using SparseArrays
using Yao: mat, ArrayReg
using Plots

# # Thermal State Properties

# We start by probing the thermal state properties of the Rydberg Hamiltonian in a 1D system. 
# Let's use the 1D chain for simplicity and vary the parameters of the Rydberg Hamiltonian and calculate the corresponding thermal state properties.
# Here, we consider a chain with 9 atoms, where nearby atoms are seperated by a distance of 5.72 μm. 
# Please refer to the [Rydberg Blockade](@ref blockade) page on tips for setting the separation distance for the atoms in preparing different ordered states.
# One can generate the system as follows using the function [`generate_sites`](@ref):

nsites = 9
atoms = generate_sites(ChainLattice(), nsites, scale = 5.72)

# We fix the Rabi frequency to be ``Ω = 2π * 4`` MHz, and study the thermal state as a function of the detuning ``Δ`` at various temperatures specified by β:

Ω = 2π * 4
Δ_step = 30
Δ = collect(LinRange(-2π * 10, 2π * 10, Δ_step));

βs = [0.00005, 0.0005, 0.005, 0.05, 0.5, 5.0, 50.0, 500.0]
labels = ["β=0.00005", "β=0.0005", "β=0.005", "β=0.05", "β=0.5", "β=5.0", "β=50.0", "β=500.0"]

all_energies = zeros(length(βs), Δ_step)
all_densities = zeros(length(βs), Δ_step, nsites)

for ii in 1:Δ_step
    h_ii = rydberg_h(atoms; Δ = Δ[ii], Ω) 
    h_m = Matrix(mat(h_ii)) 
    energies, vecs = LinearAlgebra.eigen(h_m) 
    eigenstates = [ArrayReg(vecs[:,i]) for i in eachindex(energies)] # creates the collection of eigenstates as Yao Registers

    for (β_index,β) in enumerate(βs)
        w = exp.(-β .* (energies .- energies[1]))
        all_energies[β_index, ii] = sum(w .* energies) / sum(w)

        for jj in 1:nsites
            n = rydberg_density.(eigenstates, jj) # measure the density of Rydberg excitations on each site
            all_densities[β_index, ii, jj] = sum(w .* n) / sum(w)
        end
    end 
end

# ENERGY

# Initiate plot
plt = plot(Δ / 2π, all_energies[1,:],label=labels[1])             
# uses plot() instead of plot!() as in next beta-values
xlabel!("Δ/2π (MHz) ")
ylabel!("Energy")

# Add plots
N = length(βs) - 1
for i in 1:N
    plt = plot!(Δ / 2π, all_energies[i+1,:],label=labels[i+1])
end

display(plt)
savefig("Energy_DetuningSweep.png")

# DENSITY

# We plot the densities at 3 points during the sweep: large negative detuning, zero detuning and large positive detuning.
# We compare these three plots for a high temperature and low temperature case.

# Let's start with low temperature, e.g. β = 500.0.

bar(all_densities[8,1,:], label="β=500.0")
xlabel!("Sites")
ylabel!("Rydberg density")
title!("Density Profile: 1D Chain, Δ = -2π * 10 MHz")
savefig("DensityProfile_LargeNegativeDetuning_LowT.png")

bar(all_densities[8,15,:], label="β=500.0")
xlabel!("Sites")
ylabel!("Rydberg density")
title!("Density Profile: 1D Chain, Δ = 0 MHz")
savefig("DensityProfile_ZeroDetuning_LowT.png")

bar(all_densities[8,30,:], label="β=500.0")
xlabel!("Sites")
ylabel!("Rydberg density")
title!("Density Profile: 1D Chain, Δ = 2π * 10 MHz")
savefig("DensityProfile_LargePositiveDetuning_LowT.png")


# We can nicely observe a transition into the Z2 phase throughout the sweep. This means our new code reproduces the behavior of our previous code for low temperature. (See tutorial on Adiabatic Evolution where we investigated ground state properites.)

# Now, we look a high temperature case, e.g. β = 0.00005.

bar(all_densities[1,1,:], label="β=0.00005")
xlabel!("Sites")
ylabel!("Rydberg density")
title!("Density Profile: 1D Chain, Δ = -2π * 10 MHz")
savefig("DensityProfile_LargeNegativeDetuning_HighT.png")

bar(all_densities[1,15,:], label="β=0.00005")
xlabel!("Sites")
ylabel!("Rydberg density")
title!("Density Profile: 1D Chain, Δ = 0 MHz")
savefig("DensityProfile_ZeroDetuning_HighT.png")

bar(all_densities[1,30,:], label="β=0.00005")
xlabel!("Sites")
ylabel!("Rydberg density")
title!("Density Profile: 1D Chain, Δ = 2π * 10 MHz")
savefig("DensityProfile_LargePositiveDetuning_HighT.png")

# We see that the density profile is now insensitive to the detuning.

# Insert more physics explanations.

