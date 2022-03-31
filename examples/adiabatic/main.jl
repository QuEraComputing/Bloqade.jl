# # Background

# In this example, we will show how to use the Emulator to prepare ordered ground states in Rydberg systems. 
# The example is based on experimental works in a [1D system](https://www.nature.com/articles/nature24622), and [2D system](https://www.nature.com/articles/s41586-021-03582-4). 

# Due to the strong Rydberg interactions, only one Rydberg excitation is allowed within the blockade radius [Blockade](@ref). With positive detunings, more Rydberg excitations 
# are favored (to lower the ground state(s) energy). With the interplay of these two mechanisms, different ordered states are supported depending on the strength of blockade radius and the detuning,
# such as the [``Z_N`` ordered states](https://www.nature.com/articles/nature24622) (1D) and the checkboard phase, star phase, and pure quantum phase ([straited phase](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.124.103601)) (2D). 

# We can use Quantum Annealing (QA) to prepare these quantum many-body ground states. In this process, we start with all atoms in the ground state 
# ``| 0 \rangle`` (the ground states of the many-body Hamiltonian with large negative detuning ``\Delta``). 
# Then, the Rabi frequency ``\Omega`` is turned on, and the detuning strength is ramped up from a large negative value to postive values. If this process is slow enough, the quantum state of the system stays close to the ground state of the 
# instantaneous Hamiltonian. At the end of the process, we arrive at a target Hamiltonian, and correspondingly, the prepared state is approximately the ground state for the final Hamiltonian.  


# We start by importing required libraries

using EaRyd
using CairoMakie
using EaRydPlots
using KrylovKit
using SparseArrays

# # Ground state properties

# To get a sense of the quantum phase transition, we will first vary the parameters of the Hamiltonian and calculate the corresponding ground state properties. 
# We start with a 1D chain for simplicity: we consider a chain with 9 atoms, where nearby atoms are seperated by a distance of 5.72 ``\mu m``. 
# We can generate the system as follows

nsites = 9
atoms = generate_sites(ChainLattice(), nsites, scale=5.72)

# We fix the Rabi frequency to be ``Ω = 2π * 4``, and study the ground state as a function of detuning ``Δ``. 

Ω = 2π * 4
Δ_step = 30
Δ = LinRange(-10*2π , 10*2π , Δ_step);

# We compute the Rydberg density profile for each parameter of ``\Delta``. 

density_g = zeros(Δ_step, nsites)

for ii=1: Δ_step

    h_ii = rydberg_h(atoms; Δ=Δ[ii], Ω)
    h_m = SparseMatrixCSC(h_ii)
    vals, vecs, info = KrylovKit.eigsolve(h_m,  1, :SR)
    g_state = ArrayReg(vecs[1]) 

    for jj in 1:nsites
        density_g[ii, jj] = real(expect(put(nsites, jj=>Op.n), g_state))
    end

end

# To compare, we first plot the density profile when ``\Delta= -10*2π``, 

barplot(1:nsites, density_g[1, :],axis = (xticks = 1:9,title = "Density Profile: 1D Chain, Δ = -10 * 2π", xlabel = "Sites", ylabel = "Rydberg density"))

        
# We can see that the Rydberg densities in this case is close to 0 for all sites. In contrast, for ``\Delta= 10*2π``, the density shows a clear ``Z_2`` ordered profile 
barplot(1:nsites, density_g[20, :],axis = (xticks = 1:9,title = "Density Profile: 1D Chain, Δ = 10 * 2π", xlabel = "Sites", ylabel = "Rydberg density"))
        

# More generally, we can plot an order parameter as a function of ``\Delta`` to clearly see the onset of phase transition. We can define the parameter by the difference of
# Rydberg densities in even and odd sites. 

order_para = map(1: Δ_step) do ii

    sum(density_g[ii, 1:2:8])-  sum(density_g[ii, 2:2:8])

end


lines(collect(Δ)/2π, order_para, axis = (xlabel = "Δ/2π (MHz) ", ylabel = "Order parameter"))
        
# From the density profile of ground states and the change in order parameter, we observe a phase transition with the change in ``\Delta``. 
# Below, we show that by slowly changing the parameters in the Hamiltonian, we can follow the trajectory of the ground states and adibatically evolve the atoms from the ground state to the ``Z_2`` 
# ordered state.  



# # Ordered state prepration in 1D 

# We first specify the adiabatic pulse sequence for Rabi frequency by using the built-in waveform function [`piecewise_linear`](@ref)
total_time = 3.0;
Ω_max = 2π * 4;
Ω = piecewise_linear(clocks=[0.0, 0.1, 2.1, 2.2, total_time], values=[0.0, Ω_max, Ω_max, 0, 0]);


# The detuning sequence could also be created in a similar way
U1 = -2π * 10;
U2 = 2π * 10;
Δ = piecewise_linear(clocks=[0.0, 0.6, 2.1, total_time], values=[U1, U1, U2, U2]);
    
# We plot the two waveforms:
fig = Figure();
draw!(fig[1,1], Ω, title="Ω")
draw!(fig[2,1], Δ, title="Δ")
fig

# We then generate the positions of a 1D atomic chain by using the function [`generate_sites`](@ref) 

nsites = 9
atoms = generate_sites(ChainLattice(), nsites, scale=5.72)

# Note that we have specified the nearest-neighbor atoms to be seperated by 5.72 ``\mu m`` in order to prepare a ``Z_2`` ordered state. 
# With the waveforms and atomic coordinates specified, the time-dependent Hamiltonian can be simply generated by

h = rydberg_h(atoms; Δ, Ω)

# We then specify all atoms to be in ground state initially, and set up the emulation problem by choosing the ODE solver

reg = zero_state(9);
prob = SchrodingerProblem(reg, total_time, h; dt=1e-5, adaptive=false);
integrator = init(prob, Vern8());
        
# We measure the Rydberg density for each site and time step

densities = []
for _ in integrator
    push!(densities, [expect(put(nsites, i=>Op.n), reg) for i in 1:nsites])
end
D = hcat(densities...);

# Finally, we plot the time-dependent dynamics of Rydberg density for each site

clocks = [t for t in 0:1e-5:total_time]
heatmap(clocks, 1:9, D'; axis=(xlabel="iterations", ylabel="rydberg density per site"))

        
# We can clearly see that a ``Z_2`` ordered state has been generated by the specified adiabatic pulse sequence. 
# We can also confirm this by plotting out the bitstring distribution at the final time step

bitstring_histgram(reg; nlargest=20)

# To prepare ``Z_3`` (``Z_4``) states, we can reduce the seperation between nearby atoms to 3.57 ``\mu m`` (2.87 ``\mu m``). 

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
bitstring_histgram(prob.reg; nlargest=20)


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
    
fig = Figure();
draw!(fig[1,1], Ω, title="Ω")
draw!(fig[2,1], Δ, title="Δ")
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
heatmap(clocks, 1:9, D'; axis=(xlabel="iterations", ylabel="rydberg density per site"))
