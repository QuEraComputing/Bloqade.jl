
# # Background

# The experimental study [H. Bernien, et al.](https://www.nature.com/articles/nature24622) finds that if one starts with a 
# with a particular initial state (e.g. the Neel state), the Rydberg blockade constraint results into persistent revivals of quantum dynamics.
# Later, theoretical studies (e.g. [C. J. Turner, et al.](https://www.nature.com/articles/s41567-018-0137-5)) reveal that this behavior is due to very 
# specific eigenstates embeded in the quantum many-body spectuum, called quantum many-body scars. 

# Quantum many-body scars are in anology with clasical scars in single-particle quantum chaos, where scars represent a concentration of some eigenfunctions 
# along the trajectory of classical periodic orbits. Similarly, in the quantum many-body case, the initial Neel state has a large component of these specific scar states. 
# Under the time evolution of the Rydberg Hamiltonian, the initial state undergoes the trajectory of periodic quantum orbits. The non-thermal behavior is mainly caused by such non-ergodicity 
# in Hilbert space. 

# In this example, we use the Rydberg Emulator to simulate the evolution of a fully coherent, 
# strongly interacting Rydberg system of 9 qubits.  We demonstrate the persistent revivals of many-body dynamics with measurements of the Rydberg density, 
# and entanglement entropy. For a comprehensive review of quantum many-body scars, we refer readers to this paper [M. Serbyn et al.](https://www.nature.com/articles/s41567-021-01230-2)


# We start by importing required libraries


using EaRyd
using Random
using CairoMakie


# # Build Hamiltonian

# We build a 1D-Chain with 9-atom arrangement, with each atom separated from its neighbor by 5.72 ``\mu m``. This results in a nearest-neighbor 
# interaction strength of ``2 \pi * 24`` MHz. This is much larger than the Rabi oscillations ``\Omega``, which we specify below. So the nerest-neighbor
# Rydberg atoms are within the blockade radius, such that both of the atoms can not be excited simultaneously. 

nsites = 9
atoms = generate_sites(ChainLattice(), nsites, scale=5.72)

# The waveforms are composed by two parts. For the first part, we use the adibatic evolution to prepare an ordered Neel state (see [adiabatic](@ref) for more details).

Δ1= piecewise_linear(clocks=[0.0, 0.3, 1.6, 2.2], values=[-10*2π, -10*2π, 10*2π, 10* 2π])
Ω1= piecewise_linear(clocks=[0.0, 0.05, 1.6, 2.2], values=[0.0, 4*2π, 4*2π, 0])


# The second part of the waveform has constant value of parameters, so we can use [`constant`](@ref) to construct
Ω2 = constant(duration=2.0, value=2* 2π)
Δ2 = constant(duration=2.0, value=0)

# The waveform for the whole evolution is composed by append the second part to the first part

Ω_tot = append(Ω1, Ω2)
Δ_tot = append(Δ1, Δ2)

# Note that the total evolution is 4.2 ``\mu s``.
# We then build the Hamiltonian by importing the defined lattice structure and waveforms 

h = rydberg_h(atoms; C=2 * pi * 858386, Δ=Δ_tot, Ω=Ω_tot)


# # Emulate the problem

# We evaluate the quench dynamics of the Rydberg atom array initially prepared in ground. 
# Such an initial state can be created by

reg = zero_state(9)

# We can now set up discrete time evolution problem using the ODE solver. 
prob = SchrodingerProblem(reg, 4.2, h)
integrator = init(prob, Vern8())
# Then we measure the real-time expectation value of Rydberg density, and entanglement entropy. 

entropy = Float64[]
densities = []
for _ in TimeChoiceIterator(integrator, 0.0:1e-3:4.2)
    push!(densities, [expect(put(nsites, i=>Op.n), reg) for i in 1:nsites])
    rho = density_matrix(reg, (1,2,3,4,5))
    push!(entropy, von_neumann_entropy(rho))
end


# # Plot the results 
# Now we first plot the Rydberg density for each site as a function of time

clocks = [t for t in 0:1e-3:4.2]
D = hcat(densities...)
heatmap(clocks, 1:9, D'; axis=(xlabel="iterations", ylabel="rydberg density per site"))

# We can see that the state evolve to Neel state after the first part of pulse. After that, there is a clear oscillations between the two partterns of the Rydberg density.

# We can also plot the entanglement as a function of time, 

fig = Figure(size=(5, 3));
ax = Axis(fig[1, 1])
lines!(clocks, entropy)
fig



# # A different initial state 

# In order to show that the revivals depends strongly on the initial state, 
# we now choose a different initial state, and use ['KrylovEvolution']@(ref) solver to emulate the problem  


h1= rydberg_h(atoms;C = 2π * 858386, Ω=4π)

iteration = 1:120
ts = [0.01 for _ in iteration];
hs = [h1 for _ in iteration];
clocks = cumsum(ts);
init1 = product_state(bit"100000101")
prob1 = KrylovEvolution(init1, ts, hs)
density_mat1 = zeros(nsites, length(iteration)) 

for info in prob1
    for i in 1:nsites
        density_mat1[i, info.step] = expect(put(nsites, i=>Op.n), info.reg)
    end
end

heatmap(clocks, 1:nsites, density_mat1')

# From the above figure, we see that the density does not show long-lived oscillations. 