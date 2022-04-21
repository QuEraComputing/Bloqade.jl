
# # Background

# The experimental study [H. Bernien, et al.](https://www.nature.com/articles/nature24622) finds that if one starts with a 
# with a particular initial state (e.g. the Neel state), the Rydberg blockade constraint results into persistent revivals of quantum dynamics.
# Later, theoretical studies (e.g. [C. J. Turner, et al.](https://www.nature.com/articles/s41567-018-0137-5)) reveal that this behavior is due to very 
# specific eigenstates embeded in the quantum many-body spectrum, called quantum many-body scars. 

# Quantum many-body scars are analogous to clasical scars in single-particle quantum chaos, where scars represent a concentration of some eigenfunctions 
# along the trajectory of classical periodic orbits. Similarly, in the quantum many-body case, the initial Neel state has a large component of these specific scar states. 
# Under the time evolution of the Rydberg Hamiltonian, the initial state undergoes the trajectory of periodic quantum orbits. The non-thermal behavior is mainly caused by such non-ergodicity 
# in Hilbert space. 

# In this example, we use the Rydberg Emulator to simulate the evolution of a fully coherent, 
# strongly interacting Rydberg system.  We demonstrate the persistent revivals of many-body dynamics with measurements of the Rydberg density, 
# and entanglement entropy. For a comprehensive review of quantum many-body scars, we refer readers to this paper [M. Serbyn et al.](https://www.nature.com/articles/s41567-021-01230-2)


# We start by importing required libraries

using Bloqade
using BloqadePlots
using PythonCall
using Random

plt = pyimport("matplotlib.pyplot")
# # Rabi oscillations with Rydberg blockade

# We first demonstrate that the strong Rydberg interactions have important effects on the Rabi oscillations of Rydberg atoms. 
# To do so, we consider a system with 1, 2 and 3 atoms. All the atoms are placed withint the blockade radius of any other atom (see [blockade](@ref) for more details). 

atom1 = generate_sites(ChainLattice(), 1, scale=3.0)
atom2 = generate_sites(ChainLattice(), 2, scale= 3.0)
atom3 = generate_sites(ChainLattice(), 3, scale= 3.0)

# Then a resonant Rabi driving is applied to each of the system. The Hamiltonians can be simply constructed by 
h1 = rydberg_h(atom1; Δ=0, Ω=2π*2)
h2 = rydberg_h(atom2; Δ=0, Ω=2π*2)
h3 = rydberg_h(atom3; Δ=0, Ω=2π*2)

# The initial states are chosen such that all atoms start from the ground state
reg1 = zero_state(1)
reg2 = zero_state(2)
reg3 = zero_state(3)

# We first emulate the dynamics for the single atom's case, where the intial state is quenched under a Hamiltonain with constant Rabi frequency
clocks = 0.0:1e-2:1.5
prob1 = KrylovEvolution(reg1, clocks, h1)
density_mat1 = zeros(1, length(clocks)-1)

for info in prob1
    for i in 1:1
        density_mat1[i, info.step] = expect(put(1, i=>Op.n), info.reg)
    end
end

# The Rydberg density of this atom exihibits Rabi oscillations as a function of time, shown by the plot 
fig, ax = plt.subplots()
ax.plot(clocks[1:end-1], density_mat1[1, :])
ax.set_xlabel("Time (μs)")
ax.set_ylabel("Single Rydberg Probability")
ax.set_title("Rydberg Density: Single Atom Case")
fig

# For the case of 2 and 3 atoms, if they are seperated far enough with negligible interactions, the total Rydberg excitation densities are simply the sum of each each atom. 
# However, we will show that this is not the case for systems when atoms are close to each other (which results in strong Rydberg interactions). 
# Similar to the 1 atom case, we can emulate the dynamics and get the time-dependent dynamics for each atom
prob2 = KrylovEvolution(reg2, clocks, h2);
density_mat2 = zeros(2, length(clocks)-1); 

for info in prob2
    for i in 1:2
        density_mat2[i, info.step] = expect(put(2, i=>Op.n), info.reg)
    end
end

prob3 = KrylovEvolution(reg3, clocks, h3);
density_mat3 = zeros(3, length(clocks)-1); 

for info in prob3
    for i in 1:3
        density_mat3[i, info.step] = expect(put(3, i=>Op.n), info.reg)
    end
end


# The total Rydberg density for the 1-, 2-, and 3-atom system is plotted below
fig, ax = plt.subplots()
ax.plot(clocks[1:end-1], density_mat1[1, :])
ax.set_xlabel("Time (μs)")
ax.set_ylabel("Single Rydberg Probability")
ax.set_title("Rydberg Density: Single Atom Case")
fig

# 2-atom system 
density2 = sum(density_mat2, dims=1)

fig, ax = plt.subplots()
ax.plot(clocks[1:end-1], density2[1, :])
fig

# 3-atom system 
density3 = sum(density_mat3, dims=1)

fig, ax = plt.subplots()
ax.plot(clocks[1:end-1], density3[1, :])
fig

# From the above plot, we can see that the total Rydberg density for 2 (3) atom case does not exceed 1. This is because
# it is energitically unfavorable to have more than 1 excitation due to the strong Rydberg interactions. Furthermore, the frequency of Rabi oscillation
# for the whole system depends strongly on the number of atoms. This again validates the fact that interaction plays an important role in the system's dynamics. 


# Below, we show that for a system with 9 atoms where only nearest atoms are within each other's blockade radius, the system can also exhibit nontrivial dynamics for certain initial state.  


# # Build 9-sites Hamiltonian

# We build a 1D-Chain with 9-atom arrangement, with each atom separated from its neighbor by 5.72 ``\mu m``. This results in a nearest-neighbor 
# interaction strength of ``2 \pi * 24`` MHz. This is much larger than the Rabi oscillations ``\Omega``, which we specify below. So the nerest-neighbor
# Rydberg atoms are within the blockade radius, such that both of the atoms can not be excited simultaneously. 

nsites = 9
atoms = generate_sites(ChainLattice(), nsites, scale=5.72)

# The waveforms are composed by two parts. For the first part, we use the adibatic evolution to prepare an ordered Neel state (see [adiabatic](@ref) for more details).

Δ1= piecewise_linear(clocks=[0.0, 0.3, 1.6, 2.2], values=[-10*2π, -10*2π, 10*2π, 10* 2π]);
Ω1= piecewise_linear(clocks=[0.0, 0.05, 1.6, 2.2], values=[0.0, 4*2π, 4*2π, 0]);

# The second part of the waveform has constant value of parameters, so we can use [`constant`](@ref) to construct
Ω2 = constant(duration=2.0, value=2* 2π);
Δ2 = constant(duration=2.0, value=0);

# The waveform for the whole evolution is composed by appending the second part to the first part

Ω_tot = append(Ω1, Ω2);
Δ_tot = append(Δ1, Δ2);

fig, (ax1, ax2) = plt.subplots(ncols=2)
draw!(ax1, Ω_tot)
draw!(ax2, Δ_tot)
fig

# Note that the total evolution is 4.2 ``\mu s``.
# We then build the Hamiltonian by importing the defined lattice structure and waveforms. 

h = rydberg_h(atoms; Δ=Δ_tot, Ω=Ω_tot)


# # Emulate the problem

# We evaluate the quench dynamics of the Rydberg atom array (initially prepared in the ground state). 
# The initial state can be created by

reg = zero_state(9)

# We can now set up discrete time evolution problem using the ODE solver. 
prob = SchrodingerProblem(reg, 4.2, h);
integrator = init(prob, Vern8());
# Then we measure the real-time expectation value of Rydberg density and entanglement entropy. 

entropy = Float64[]
densities = []
for _ in TimeChoiceIterator(integrator, 0.0:1e-3:4.2)
    push!(densities, [expect(put(nsites, i=>Op.n), reg) for i in 1:nsites])
    rho = density_matrix(reg, (1,2,3,4,5))
    push!(entropy, von_neumann_entropy(rho))
end


# # Plot the results 
# We first plot the Rydberg density for each site as a function of time

clocks = [t for t in 0:1e-3:4.2]
D = hcat(densities...)

fig, ax = plt.subplots(figsize = (10,4))
ax.imshow(real(D), interpolation="nearest", aspect="auto")
ax.set_xlabel("iterations")
ax.set_ylabel("rydberg density per site")
fig

# We can see that the state evolves to a Neel state after the first part of pulse. After that, there are clear oscillations between the two partterns of the Rydberg density.

# We can also plot the entanglement as a function of time

fig, ax = plt.subplots(figsize = (10,4))
ax.plot(clocks, entropy)
ax.set_xlabel("Time (μs)")
ax.set_ylabel("Entanglement Entropy")
fig

# # A different initial state 

# In order to show that the revivals depends strongly on the initial state, 
# we now choose a different initial state, and use the ['KrylovEvolution']@(ref) solver to emulate the problem  


hd = rydberg_h(atoms;Ω=4π)
clocks = 0.0:1e-2:1.2

init_d = product_state(bit"100000101")
prob_d = KrylovEvolution(init_d, clocks, hd)
density_mat_d = zeros(nsites, length(clocks)-1) 

for info in prob_d
    for i in 1:nsites
        density_mat_d[i, info.step] = expect(put(nsites, i=>Op.n), info.reg)
    end
end

fig, ax = plt.subplots(figsize = (10,4))
ax.imshow(real(density_mat_d), interpolation="nearest", aspect="auto")
ax.set_xlabel("iterations")
ax.set_ylabel("rydberg density per site")
fig
        
# From the above figure, we see that the density does not show long-lived oscillations. 
