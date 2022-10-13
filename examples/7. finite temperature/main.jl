# # Adiabatic Evolution
# ## Background

# In this example, we will show how to use Bloqade to prepare ordered ground states in the Rydberg system. 
# The example is based on the experimental works in a [1D system](https://www.nature.com/articles/nature24622) and [2D system](https://www.nature.com/articles/s41586-021-03582-4).
# The Rydberg Hamiltonian can be found in the [Bloqade](@ref) page.

# Due to the strong Rydberg interactions, only one Rydberg excitation is allowed within the blockade radius (see [Rydberg Blockade](@ref blockade)). With a positive detuning Δ, more Rydberg excitations 
# are favored (to lower the ground state(s) energy). The interplay of these two mechanisms allows the creation of different ordered states depending on the strength of the blockade radius and the detunings,
# such as the [``Z_N`` ordered states](https://www.nature.com/articles/nature24622) in 1D and the checkerboard phase, the star phase, and a pure quantum phase (the striated phase) in 2D 
# (see the [experimental](https://www.nature.com/articles/s41586-021-03582-4) and [theory](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.124.103601) papers).

# Here, we use the Quantum Adiabatic Algorithm (QAA) to prepare these quantum many-body ground states. To do that, we can start with all atoms in the ground state 
# ``| 0 \rangle``, which is the ground state of the many-body Hamiltonian with a large negative detuning ``\Delta``. 
# Then, the Rabi frequency ``\Omega`` is turned on, and the detuning strength is ramped up from a large negative value to postive values. If this process is slow enough, the quantum state of the system stays close to the ground state of the 
# instantaneous Hamiltonian. At the end of this process, we arrive at a target Hamiltonian, and correspondingly, the prepared state is approximately the ground state for the final Hamiltonian.
# A quantum phase transition typically occurs during this process and one can probe the phase transition and critical phenomena by simulating and understanding the quantum dynamics.

# Let's start by importing the required libraries:

using Bloqade
using PythonCall
using LinearAlgebra
using SparseArrays

plt = pyimport("matplotlib.pyplot");

# # Ground State Properties

# We start by probing the ground state properties of the Rydberg Hamiltonian in a 1D system. 
# Let's use the 1D chain for simplicity and vary the parameters of the Rydberg Hamiltonian and calculate the corresponding ground state properties.
# Here, we consider a chain with 9 atoms, where nearby atoms are seperated by a distance of 5.72 μm. 
# Please refer to the [Rydberg Blockade](@ref blockade) page on tips for setting the separation distance for the atoms in preparing different ordered states.
# One can generate the system as follows using the function [`generate_sites`](@ref):

nsites = 9
atoms = generate_sites(ChainLattice(), nsites, scale = 5.72)

# We fix the Rabi frequency to be ``Ω = 2π * 4`` MHz, and study the ground state as a function of the detuning ``Δ``:

β = 0.05
Ω = 2π * 4
Δ_step = 30
Δ = collect(LinRange(-2π * 10, 2π * 10, Δ_step));

# The Rydberg density profile can be computed for each parameter of ``\Delta`` as: 

density_g = zeros(Δ_step, nsites)

for ii in 1:Δ_step
    h_ii = rydberg_h(atoms; Δ = Δ[ii], Ω) # create the Rydberg Hamiltonian
    h_m = Matrix(mat(h_ii)) # convert the Hamiltonian into a matrix
    energies, vecs = LinearAlgebra.eigen(h_m) # find the full spectrum
    eigenstates = [ArrayReg(vecs[:,i]) for i in eachindex(energies)] # creates the collection of eigenstates as Yao Registers
    w = exp.(-β .* (energies .- energies[1]))

    
    for jj in 1:nsites
        n = rydberg_density.(eigenstates, jj) # measure the density of Rydberg excitations on each site
        density_g[ii, jj] = sum(w .* n) / sum(w)
    end
end

# To compare, we first plot the density profile when ``\Delta= -2π * 10`` MHz: 

fig, ax = plt.subplots(figsize = (10, 4))
ax.bar(1:nsites, density_g[1, :])
ax.set_xticks(1:nsites)
ax.set_xlabel("Sites")
ax.set_ylabel("Rydberg density")
ax.set_title("Density Profile: 1D Chain, Δ = -2π * 10 MHz")
fig

# We can see that the Rydberg densities in this case is close to 0 for all sites. In contrast, for ``\Delta= 2π * 10`` MHz, the density shows a clear ``Z_2`` ordered profile:
fig, ax = plt.subplots(figsize = (10, 4))
ax.bar(1:nsites, density_g[30, :])
ax.set_xticks(1:nsites)
ax.set_xlabel("Sites")
ax.set_ylabel("Rydberg density")
ax.set_title("Density Profile: 1D Chain, Δ = 2π * 10 MHz")
fig

# More generally, we can plot an order parameter as a function of ``\Delta`` to clearly see the onset of phase transition. 
# The order parameter can be defined as the difference of Rydberg densities on even and odd sites. 

order_para = map(1:Δ_step) do ii
    return sum(density_g[ii, 1:2:nsites]) - sum(density_g[ii, 2:2:nsites])
end

fig, ax = plt.subplots(figsize = (10, 4))
ax.plot(Δ / 2π, order_para)
ax.set_xlabel("Δ/2π (MHz) ")
ax.set_ylabel("Order parameter")
fig

# From the density profile of ground states and the change in the order parameter, we can observe a phase transition with changing ``\Delta``. 
# Below, we show that by slowly changing the parameters of the Hamiltonian, we can follow the trajectory of the ground states and adiabatically evolve the atoms from the ground state to the ``Z_2`` 
# ordered state.


fig, ax = plt.subplots(figsize = (10, 4))
shw = ax.imshow(transpose(real(density_g)), interpolation = "nearest", aspect = "auto", extent = [0, Δ[end], 0.5, nsites + 0.5],vmin=0,vmax=1)
ax.set_xlabel("Δ/2π (MHz) ")
ax.set_ylabel("site")
ax.set_xticks(Δ[1:end:3])
ax.set_yticks(1:nsites)
bar = fig.colorbar(shw)
fig

