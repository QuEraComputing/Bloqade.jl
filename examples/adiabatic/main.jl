# # Background

# In this example, we will show how to use Bloqade to prepare ordered ground states in the Rydberg system. 
# The example is based on the experimental works in a [1D system](https://www.nature.com/articles/nature24622) and [2D system](https://www.nature.com/articles/s41586-021-03582-4).
# The Rydberg Hamiltonian can be found in [Bloqade](@ref).

# Due to the strong Rydberg interactions, only one Rydberg excitation is allowed within the blockade radius (see [Rydberg Blockade](@ref)). With a positive detuning Δ, more Rydberg excitations 
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
using BloqadePlots
using KrylovKit
using SparseArrays

plt = pyimport("matplotlib.pyplot")

# # Ground state properties

# We start by probing the ground state properties of the Rydberg Hamiltonian in a 1D system. 
# Let's use the 1D chain for simplicity and vary the parameters of the Rydberg Hamiltonian and calculate the corresponding ground state properties.
# Here, we consider a chain with 9 atoms, where nearby atoms are seperated by a distance of 5.72 ``\mu m``. 
# One can generate the system as follows:

nsites = 9
atoms = generate_sites(ChainLattice(), nsites, scale = 5.72)

# We fix the Rabi frequency to be ``Ω = 2π * 4`` MHz, and study the ground state as a function of the detuning ``Δ``.

Ω = 2π * 4
Δ_step = 30
Δ = LinRange(-2π * 10, 2π * 10, Δ_step);

# The Rydberg density profile can be computed for each parameter of ``\Delta``. 

density_g = zeros(Δ_step, nsites)

for ii in 1:Δ_step
    h_ii = rydberg_h(atoms; Δ=Δ[ii], Ω) # create the Rydberg Hamiltonian
    h_m = mat(h_ii) # convert the Hamiltonian into a matrix
    vals, vecs, info = KrylovKit.eigsolve(h_m,  1, :SR) # find the ground state eigenvalue and eigenvector
    g_state = ArrayReg(vecs[1]) # creates the initial state with all atoms in ``| 0 \rangle`` state

    for jj in 1:nsites
        density_g[ii, jj] = real(expect(put(nsites, jj=>Op.n), g_state)) # measure the density of Rydberg excitations on each site
    end
end

# To compare, we first plot the density profile when ``\Delta= -2π * 10`` MHz, 

fig, ax = plt.subplots(figsize = (10, 4))
ax.bar(1:nsites, density_g[1, :])
ax.set_xticks(1:nsites)
ax.set_xlabel("Sites")
ax.set_ylabel("Rydberg density")
ax.set_title("Density Profile: 1D Chain, Δ = -2π * 10")
fig

# We can see that the Rydberg densities in this case is close to 0 for all sites. In contrast, for ``\Delta= 2π * 10`` MHz, the density shows a clear ``Z_2`` ordered profile
fig, ax = plt.subplots(figsize = (10, 4))
ax.bar(1:nsites, density_g[30, :])
ax.set_xticks(1:nsites)
ax.set_xlabel("Sites")
ax.set_ylabel("Rydberg density")
ax.set_title("Density Profile: 1D Chain, Δ = -2π * 10 MHz")
fig
        

# More generally, we can plot an order parameter as a function of ``\Delta`` to clearly see the onset of phase transition. 
# The order parameter can be defined as the difference of Rydberg densities on even and odd sites. 

order_para = map(1: Δ_step) do ii
    sum(density_g[ii, 1:2:nsites]) - sum(density_g[ii, 2:2:nsites])
end

fig, ax = plt.subplots(figsize = (10,4))
ax.plot(Δ/2π, order_para)
ax.set_xlabel("Δ/2π (MHz) ")
ax.set_ylabel("Order parameter")
fig
        
# From the density profile of ground states and the change in the order parameter, we can observe a phase transition with changing ``\Delta``. 
# Below, we show that by slowly changing the parameters of the Hamiltonian, we can follow the trajectory of the ground states and adiabatically evolve the atoms from the ground state to the ``Z_2`` 
# ordered state.


# # Preparation of ordered states in 1D 

# We first specify the adiabatic pulse sequence for Rabi frequency by using the built-in waveform function [`piecewise_linear`](@ref).
total_time = 3.0;
Ω_max = 2π * 4;
Ω = piecewise_linear(clocks = [0.0, 0.1, 2.1, 2.2, total_time], values = [0.0, Ω_max, Ω_max, 0, 0]);


# The detuning sequence can also be created in a similar way.
U1 = -2π * 10;
U2 = 2π * 10;
Δ = piecewise_linear(clocks = [0.0, 0.6, 2.1, total_time], values = [U1, U1, U2, U2]);
    
# We plot the two waveforms:
fig, (ax1, ax2) = plt.subplots(ncols = 2, figsize = (12, 4))
draw!(ax1, 1/2π .* Ω)
ax1.set_xlabel("time (μs)")
ax1.set_ylabel("Ω/2π MHz")
draw!(ax2, 1/2π .* Δ)
ax2.set_xlabel("time (μs)")
ax2.set_ylabel("Δ/2π MHz")
fig

# We then generate the positions of a 1D atomic chain by using the function [`generate_sites`](@ref) 

nsites = 9
atoms = generate_sites(ChainLattice(), nsites, scale=5.72)

# Note that we have specified the nearest-neighbor atoms to be seperated by 5.72 ``\mu m`` in order to prepare a ``Z_2`` ordered state. 
# With the waveforms and atomic coordinates specified, the time-dependent Hamiltonian can be simply generated by

h = rydberg_h(atoms; Δ, Ω)

# We then specify all atoms to be in ground state initially, and set up the emulation problem by choosing the ODE solver

reg = zero_state(9);
prob = SchrodingerProblem(reg, total_time, h);
integrator = init(prob, Vern8());
        
# We measure the Rydberg density for each site and time step

densities = []
for _ in TimeChoiceIterator(integrator, 0.0:1e-3:total_time)
    push!(densities, [expect(put(nsites, i=>Op.n), reg) for i in 1:nsites])
end
D = hcat(densities...);

# Finally, we plot the time-dependent dynamics of Rydberg density for each site

clocks = [t for t in 0:1e-3:total_time]
fig, ax = plt.subplots(figsize = (10,4))
ax.imshow(real(D), interpolation="nearest", aspect="auto")
ax.set_xlabel("iterations")
ax.set_ylabel("rydberg density per site")
fig
        
# We can clearly see that a ``Z_2`` ordered state has been generated by the specified adiabatic pulse sequence. 
# We can also confirm this by plotting out the bitstring distribution at the final time step

bitstring_hist(reg; nlargest=20)

# To prepare ``Z_3`` or ``Z_4`` states, we can reduce the seperation between nearby atoms to 3.57 ``\mu m`` or 2.87 ``\mu m`` respectively. 

# # Run in the blockade subspace  

# In the above example, we have run the emulation in fullspace, which can be slow and only works for small system. We can turn the above simulation into a blockade subspace simulation 
# by changing the register to a RydbergReg by feeding a subspace object.

# The subspace can be found by looking up the independent set of the graph constructed by a blockade radius; here we choose the radius to be 6.2 ``\mu m``

space = blockade_subspace(atoms, 6.2);

# Then we create our register in subspace

reg = zero_state(space)

# The rest of code will be the same as the full space 

prob = SchrodingerProblem(reg, total_time, h)
emulate!(prob)
bitstring_hist(prob.reg; nlargest=20)


# # State preparation in 2D

#  Now we will show how to prepare a 2D checkboard phase. Most of codes will be the same as the 1D case, except that we will choose slightly different 
# parameters and specify a square lattice instead of a chain

nx, ny = 3, 3
nsites = nx*ny
atoms = generate_sites(SquareLattice(), nx, ny, scale=6.7)

total_time = 2.9
Ω_max = 2π * 4.3
Ω = piecewise_linear(clocks=[0.0, 0.3, 2.6, total_time], values=[0.0, Ω_max , Ω_max , 0]);

U = 2π * 15.0
Δ = piecewise_linear(clocks=[0.0, 0.3, 2.6, total_time], values=[-U, -U, U , U]);

fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(10, 4))
draw!(ax1, Ω)
draw!(ax2, Δ)
fig

h = rydberg_h(atoms; Δ, Ω)


reg = zero_state(9)
prob = SchrodingerProblem(reg, total_time, h)
integrator = init(prob, Vern8())
densities = [];

# If you are using adaptive steps (by default),
# you can use `TimeChoiceIterator` to specify the
# points you would like to stop by
for _ in TimeChoiceIterator(integrator, 0.0:1e-3:total_time)
    push!(densities, [expect(put(nsites, i=>Op.n), reg) for i in 1:nsites])
end
D = hcat(densities...)

clocks = [t for t in 0:1e-3:total_time]

fig, ax = plt.subplots(figsize = (10,4))
ax.imshow(real(D), interpolation="nearest", aspect="auto")
ax.set_xlabel("iterations")
ax.set_ylabel("rydberg density per site")
fig
