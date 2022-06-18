# # Simulation of lattice gauge theory with Rydberg atoms
# ## Introduction
# In previous examples, we have shown how to prepared $Z_2$ ordered ground state for the Rydberg system, and the quantum scar phenomenon, which is the oscillation between two ordered patterns of Rydberg densities. We note that these are achieved by tuning the detuning and Rabi frequency of the lasers that address all the atoms simulatenously. 

# In this tutorial, we shall simulate the dynamics of latticee gauge theory (LGT) with a 1D Rydberg atom chain. It turns out that the $Z_2$ ground state of the Rydberg chain corresponds to the vacuum state of the studied LGT, and the quantum scar is nothing but the string-inversion mechanism in the gauge theory. More interestingly, by locally detuning certain atoms, we can creat defects in the chain and simulate the propagation of particle-antiparticle pairs. This tutorial is inspired by the paper https://journals.aps.org/prx/pdf/10.1103/PhysRevX.10.021041


######## The prerequisite of simulating LGT is to prepare the so-called vacuum state for the LGT. The mapping between the Rydberg ground state and the ground state of the interested LGT is illustrated below. To put it in words, Rydberg states in the odd sites and even states correspond to vacuum $\emptyset$ and a particle $q$ respectively, whereas ground states in the odd sites and even states correspond to an anti-particle $\bar{q}$ and vacuum $\emptyset$ respectively. Hence the vacuum state of the LGT is a configuration where all the odd sites are occupied by Rydberg states whereas the even sites are ground states, which can be simulated on a Rydberg chain provided that the neighboring atoms are within the blockade radius.

######We know how to prepare the 1D AF state from the 0th use case, and starting from the ground state, we simulate the dynamics of the LGT using the Rydberg system.

# We first import the required packages 

using Bloqade
using PythonCall
using StatsBase
using Distributed
using BitBasis
using Plots
plt = pyimport("matplotlib.pyplot");


# # Mapping between the Rydberg system and the LGT

# In this tutorial, we are interested in the so-called quantum link model (QLM) formulatino of the LGT. In this formalism, depending on the configurations of the even and odd sites, the bonds between them could be interpretted as a particle $q$, an antipartile $\bar{q}$ or a vacuum $\emptyset$. More specifically, for the bond between an odd and an even sites, it corresponds to an antiparticle if both atoms are in the ground states, otherwise it is interpretted as a vacuum state; for the bond between an even and an odd sites, on the other hand, it corresponds to a vacuum if both atoms are in the ground states, otherwise it is interpretted as an antiparticle. Further, the Rydberg states at the odd (even) states are interpretted as electric fields pointing to the right (left), whereas the ground states at the odd (even) states are electric field pointing to the left (right). The mapping is summarized in the figure below (source: [Federica M. Surace, et al.](https://journals.aps.org/prx/pdf/10.1103/PhysRevX.10.021041).  


# ![mapping](../../../assets/LGT_mapping.png)


# The LGT dynamics starts from the "anti-string" state with all electric fields pointing to the right. This is nothing but the $Z_2$ ordered state in the language of Rydberg system, and we have learnt how to prepare it in previous tutorials. Here, we are interested in a 1D lattice with 21 atoms. The neighboring atoms are separated by $5.5μm$ such that they are blockaded throughout the dynamics.  

a = 5.5;
N = 21;
atoms = generate_sites(ChainLattice(), N, scale=a);
subspace = blockade_subspace(atoms, a); 

# In order to prepare the anti-string state of the LGT, we use piecewise linear waveforms for both detuning and the Rabi frequency. The waveforms will last for $3.5μs$. 

total_time = 3.5; 
Ωmax = 2π * 5;
Δmin = -2π * 10;
Δmax = 2π * 10;

Δ1 = piecewise_linear(clocks = [0.0, 0.2, total_time], values = [Δmin, Δmin, Δmax]); 
Ω1 = piecewise_linear(clocks = [0.0, 0.2, total_time-0.0001, total_time], values = [0.0, Ωmax, Ωmax, 0]);

# The waveforms for the detuning and the Rabi frequency are shown below
fig1, (ax1, ax2) = plt.subplots(ncols = 2, figsize = (12, 4))
Bloqade.plot!(ax1, Ω1)
Bloqade.plot!(ax2, Δ1)
ax1.set_ylabel("Ω1/2π (MHz)")
ax2.set_ylabel("Δ1/2π (MHz)")
ax1.grid()
ax2.grid()
fig1


# For later purpose, we define the function `Get_Ryd_density` which simulates the dynamics for a given set of detuning and Rabi frequency $(Δ, Ω)$, and returns the final state and Rydberg density. 

function Get_Ryd_density(Δ, Ω; dt=1e-3)
    h = rydberg_h(atoms; Δ=Δ, Ω=Ω)
    reg = zero_state(subspace)

    duration = Δ.duration
    prob = SchrodingerProblem(reg, duration, h, progress=true);
    integrator = init(prob, Vern8());
    densities = []
    for _ in TimeChoiceIterator(integrator, 0.0:dt:duration)
        normalize!(reg) 
        push!(densities, rydberg_density(reg)) 
    end

    return reg, densities
end

# We can confirm that the waveforms produce the desired anti-string state of the LGT, by simulating the dynamics governed by the LGT Hamiltonian, followed by plotting the density profile, as shown below

reg1, dens1 = Get_Ryd_density(Δ1, Ω1) ;
fig2, ax = plt.subplots(figsize = (10, 4)) ;
ax.bar(1:N, dens1[end]) ;
fig2 

# Recall the mapping between the Rydberg chain and the LGT illustrated above, we see that the final state is an approximation of the anti-string state of the LGT, or a $Z_2$ ordered state. It is not an ideal one where the discrepancy is more promonient at the center compared to the edge of the chain. But as we will see later, the prepared state is sufficient for demonstrating the dynamics of the interested LGT. 

# # Propagation of particle-antiparticle pairs

# Up to now, what we have done is simply reproducing the results from previous tutorials. Next, we will prepare defects in the anti-string state, which are links with right-pointing electric fields. The domain walls between anti-string and string states will host particles, whereas those between string and anti-string states will host anti-particles. These can be seen via the mapping between Rydberg system and LGT illustrated above. Interestingly, the particle and antiparticle always come in pairs, and their time evolution exhibits a light cone, in which the string-antistring oscillation is out-of-phase compared to that outside of the light cone. This is illustrated below (source: [Federica M. Surace, et al.](https://journals.aps.org/prx/pdf/10.1103/PhysRevX.10.021041)).


# ![propagation](../../../assets/LGT_particle_antiparticle_propagation.png)

# ## Site-dependent waveforms

# To realize the defects, we will apply a $\pi$-pulse to the target atoms for transitioning them from the Rydberg state to the ground state. 

Ωq = Ωmax; 
tq = pi/Ωq; 

# The waveform for the Rabi frequency for creating one and two defects are defined as following

Ω2_single_defect = map(1:length(atoms)) do idx
    if idx == floor(Int, N/2)+1
        append(Ω1, constant(duration=tq, value=Ωq))
    else
        append(Ω1, constant(duration=tq, value=0))
    end
end

Ω2_two_defects = map(1:length(atoms)) do idx ; 
    if idx == floor(Int, N/3) || idx == floor(Int, N-N/3)+1
        append(Ω1, constant(duration=tq, value=Ωq))
    else
        append(Ω1, constant(duration=tq, value=0))
    end
end

# We append to a constant waveform with zero amplitude to the detuning  such that it has the same duration as the Rabi frequencies.

Δ2 = append(Δ1, constant(duration=tq, value=0)); 

# As an example, for the case with a single defect, we show the Rabi frequency for the central site, which is the defect, and other sites separately below. 

fig3, (ax1, ax2) = plt.subplots(ncols = 2, figsize = (12, 4))
for idx in 1 : length(atoms)
    if idx == floor(Int, N/2)+1
        Bloqade.plot!(ax1, Ω2_single_defect[idx])
    else
        Bloqade.plot!(ax2, Ω2_single_defect[idx])
    end
end
ax1.grid()
ax2.grid()
ax1.set_title("Rabi frequency for the central site")
ax2.set_title("Rabi frequency for other sites")

fig3

# We can confirm that the waveforms produce the desired domain walls for the LGT states, by simulating the dynamics governed by the LGT Hamiltonian, followed by plotting the density profile. 

reg2, dens2 = Get_Ryd_density(Δ2, Ω2_single_defect)
reg3, dens3 = Get_Ryd_density(Δ2, Ω2_two_defects)

fig4, (ax1, ax2) = plt.subplots(nrows = 2, figsize = (10, 4))
ax1.bar(1:N, dens2[end])
ax2.bar(1:N, dens3[end])
fig4


# Again, we see that the Rydberg density at the defects are not exactly zero, but the prepared state, as we shall below, serve as a good initial state to study the Propagation of particle-antiparticle pairs in LGT. 

# We define the very last piece in the Rabi frequency and detuning that govern the time evolution of the Rydberg chain with defects. 
Ωq2 = Ωmax ; 
tq2 = 40/Ωq2 ; 

Ω3_single_defect = map(1:length(atoms)) do idx
    append(Ω2_single_defect[idx], constant(duration=tq2, value=Ωq2))
end       
Ω3_two_defects = map(1:length(atoms)) do idx
    append(Ω2_two_defects[idx], constant(duration=tq2, value=Ωq2))
end       

Δ3 = append(Δ2, constant(duration=tq2, value=-π));


# Again, as an example, for the case with a single defect, we show the Rabi frequency for the central site, which is the defect, and other sites separately below. 

fig5, (ax1, ax2) = plt.subplots(ncols = 2, figsize = (12, 4))
for idx in 1 : length(atoms)
    if idx == floor(Int, N/2)+1
        Bloqade.plot!(ax1, Ω3_single_defect[idx])
    else
        Bloqade.plot!(ax2, Ω3_single_defect[idx])
    end
end
ax1.grid()
ax2.grid()
ax1.set_title("Rabi frequency for the central site")
ax2.set_title("Rabi frequency for other sites")

fig5

# ## Simulation particle-antiparticle pairs in LGT dynamics 

# With the waveforms defined, we can run the simulation to evolve the Rydberg chains with defects and collect the final Rydberg densities.
densities_single_defect = Get_Ryd_density(Δ3, Ω3_single_defect)[2];
densities_two_defects = Get_Ryd_density(Δ3, Ω3_two_defects)[2];

D_single_defect = hcat(densities_single_defect...);
D_two_defects = hcat(densities_two_defects...);

# To better visualize the propagation of particle-antiparticle pairs, we shall only show the Rydberg densities starting from time when the ground state of the defect chain is prepared. 

ind0 = 3550; # 3500
D_single_defect = D_single_defect[:, ind0:end];
D_two_defects = D_two_defects[:, ind0:end];

clocks = 0:1e-3:Δ3.duration;
clocks = clocks[ind0: end];

# Then we plot the Rydberg density as a function of time, where the two panels correspond to the cases with single and two defects respectively
yticks = range(clocks[1], stop=clocks[end], length=10);
yticks = [string(ytick)[1:3] for ytick in yticks];

h1 = heatmap(real(transpose(D_single_defect)), legend=:none, xlabel="sites", ylabel="time (μs)", yticks=(range(1, length(clocks), length=10), yticks));
h2 = heatmap(real(transpose(D_two_defects)), legend=:none, xlabel="sites", yticks=(range(1, length(clocks), length=10), yticks));
l = @layout[grid(1,2) a{0.05w}]
p = plot([h1, h2]..., heatmap((0:0.01:1).*ones(101,1), legend=:none, xticks=:none, yticks=(1:10:101, string.(0:0.1:1))), layout=l)

# From the left panel, we can observe a light-cone-shaped region originating from the particle-antiparticle pair in the vacuum. At the right panel, we show the interference of two light cones, which produces an additional change of periodicity corresponding to the elastic scattering. When the particle or antiparticle reaches the boundary of the chain, it will be scattered back as observed. For more details, the interested readers are referred to the paper https://journals.aps.org/prx/pdf/10.1103/PhysRevX.10.021041. 

# In summary, we have shown that the ground state and the dynamics of certain LGT can be simulated by a 1D chain of Rydberg atoms. More interestingly, defects can be introduced by locally addressing certain atoms in the chain, and with that we can simulate the propagation of particle-antiparticle pairs in the LGT dynamics.


