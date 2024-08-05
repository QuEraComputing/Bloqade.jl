# # Simulation of a Lattice Gauge Theory with Rydberg Atoms
# ## Introduction
# In previous examples, we have shown how to prepare $Z_2$ ordered ground state for the Rydberg system, and discussed the quantum scar phenomenon. 
# We note that these are achieved by tuning the detuning and the Rabi frequency of the lasers that address all the atoms simultaneously. 

# In this tutorial, we shall simulate the dynamics of a Lattice Gauge Theory (LGT) with a 1D Rydberg atom chain. 
# In the context of gauge theories, it turns out that the $Z_2$ ground state and the quantum scar of the Rydberg chain correspond to 
# the "string" state and the string-inversion mechanism of the studied LGT respectively. More interestingly, by 
# locally addressing certain atoms, we can create defects in the chain and simulate the propagation of particle-antiparticle pairs. 
# This tutorial is inspired by the paper [F. M. Surace et al. (10.1103/PhysRevX.10.021041)](https://journals.aps.org/prx/pdf/10.1103/PhysRevX.10.021041). 


# We first import the required packages 

using Bloqade
using Bloqade.CairoMakie

# ## Mapping between the Rydberg system and the LGT

# In this tutorial, we are interested in the so-called quantum link model (QLM) formulation
# of the LGT. In this formalism, depending on the configurations of the even and odd sites,
# the bonds between them could be interpreted as a particle $q$, an antiparticle $\bar{q}$
# or a vacuum state. More specifically, the bond between an odd and an even sites corresponds
# to an antiparticle if both atoms are in the ground states, otherwise it is interpreted as
# a vacuum state. On the other hand, the bond between an even and an odd sites corresponds
# to a particle if both atoms are in the ground states, otherwise it is interpreted as a vacuum.
# Further, the Rydberg states at the odd (even) sites are interpreted as electric fields
# pointing to the left (right), whereas the ground states at the odd (even) sites are electric
# field pointing to the right (left). The electric fields correspond to the red and blue arrows
# in the following figure, which summarizes the mappings described above
# (source: [F. M. Surace et al. (10.1103/PhysRevX.10.021041)](https://journals.aps.org/prx/pdf/10.1103/PhysRevX.10.021041)).

# ![mapping](../../../assets/LGT_mapping.png)

# ## Preparing the initial state for the LGT dynamics

# The LGT dynamics starts from the "anti-string" state with all electric fields pointing
# to the right. This is nothing but the $Z_2$ ordered state in the language of Rydberg
# system, and we have seen how to prepare it in previous tutorials. Here, we are interested
# in a 1D lattice with 21 atoms. The neighboring atoms are separated by 5.5 μm such that
# they are blockaded throughout the dynamics (the typical value of the Rabi frequency is
# ``2\pi \times 5`` MHz throughout the dynamics, which corresponds to blockade radius ``R_b\approx7.46`` μm.):

a = 5.5;
N = 21;
atoms = generate_sites(ChainLattice(), N, scale=a);
subspace = blockade_subspace(atoms, a); 

# In order to prepare the anti-string state of the LGT, we use piecewise linear waveforms for both the detuning and the Rabi frequency. The waveforms will last for 3.5 μs:

total_time = 3.5; 
Ωmax = 2π * 5;
Δmin = -2π * 10;
Δmax = 2π * 10;

Δ1 = piecewise_linear(clocks = [0.0, 0.2, total_time], values = [Δmin, Δmin, Δmax]); 
Ω1 = piecewise_linear(clocks = [0.0, 0.2, total_time], values = [0.0, Ωmax, Ωmax]);

# The waveforms for the detuning and the Rabi frequency are shown below
fig1 = Figure(size=(960, 320))
ax1 = Axis(fig1[1, 1], ylabel="Ω1/2π (MHz)", xgridvisible=true, ygridvisible=true)
ax2 = Axis(fig1[1, 2], ylabel="Δ1/2π (MHz)", xgridvisible=true, ygridvisible=true)
Bloqade.plot!(ax1, Ω1)
Bloqade.plot!(ax2, Δ1)
fig1

# In order to simulate the gauge theory dynamics, we define the function
# `get_average_rydberg_densities` which takes in a detuning
# and Rabi frequency $(Δ, Ω)$, and returns the final state and the Rydberg density: 

function get_average_rydberg_densities(Δ, Ω; dt=1e-3)
    h = rydberg_h(atoms; Δ=Δ, Ω=Ω)
    reg = zero_state(subspace)

    duration = Ω.duration
    prob = SchrodingerProblem(reg, duration, h, progress=true);
    integrator = init(prob, Vern8());
    densities = []
    for _ in TimeChoiceIterator(integrator, 0.0:dt:duration)
        normalize!(reg) 
        push!(densities, rydberg_density(reg)) 
    end

    return densities
end;

# We can confirm that the waveforms produce the desired anti-string state of the LGT, 
# by simulating the dynamics governed by the waveforms, followed by plotting the density profile, as shown below:

dens1 = get_average_rydberg_densities(Δ1, Ω1) ;

fig2 = Figure(size = (800, 320))
ax = Axis(fig2[1, 1], xlabel="site index", ylabel="Average Rydberg densities")
barplot!(ax, 1:N, dens1[end])
fig2

# Recalling the mapping between the Rydberg chain and the LGT illustrated above, 
# we see that the final state is an approximation of the anti-string state of the LGT, or a $Z_2$ ordered state of the Rydberg chain. 
# It is not perfectly $Z_2$ ordered where the discrepancy is more prominent at the center compared to the edge of the chain. 
# But as we will see later, the prepared state is sufficient for demonstrating the dynamics of the interested LGT. 

# ## Propagation of Particle-Antiparticle Pairs

# Next, we prepare defects in the anti-string state, which are links with right-pointing electric fields. 
# The domain walls between anti-string and string states will host particles, whereas those between string and anti-string states will host anti-particles. 
# These can be seen via the mapping between Rydberg a system and LGT illustrated above. 
# Interestingly, the particle and antiparticle always come in pairs, and their time evolution exhibits light cones, 
# in which the string-antistring oscillation is out-of-phase compared to that outside of the light cone. 
# This is illustrated below (source: [F. M. Surace et al. (10.1103/PhysRevX.10.021041)](https://journals.aps.org/prx/pdf/10.1103/PhysRevX.10.021041)):


# ![propagation](../../../assets/LGT_particle_antiparticle_propagation.png)

# ### Site-Dependent Waveforms

# To realize the defects, we turn off the detuning for the target atoms while maintaining the same Rabi frequency for all the atoms. 
# This effectively applies a $\pi$-pulse to the target atoms for transitioning them from the Rydberg state to the ground state: 

Δq = 0.0; 
tq = π/Ωmax; 

# The waveforms of the detunings for creating one and two defects are defined as following:

Δ2_single_defect = map(1:length(atoms)) do idx
    if idx == floor(Int, N/2)+1
        append(Δ1, constant(duration=tq, value=Δq))
    else
        append(Δ1, constant(duration=tq, value=Δmax))
    end
end ; 

Δ2_two_defects = map(1:length(atoms)) do idx
    if idx == floor(Int, N/3) || idx == floor(Int, N-N/3)+1
        append(Δ1, constant(duration=tq, value=Δq))
    else
        append(Δ1, constant(duration=tq, value=Δmax))
    end
end ; 

# We append a constant waveform with the same amplitude to the Rabi frequency such that it has the same duration as the detunings:

Ω2 = append(Ω1, constant(duration=tq, value=Ωmax)); 

# As an example, for the case with a single defect, we show the detuning for the central site, which is the defect, and those for other sites separately below:

fig3 = Figure(size=(960, 320))
ax1 = Axis(fig3[1, 1], title="Detuning for the central site", xgridvisible=true, ygridvisible=true)
ax2 = Axis(fig3[1, 2], title="Detunings for other sites", xgridvisible=true, ygridvisible=true)

for idx in 1 : length(atoms)
    if idx == floor(Int, N/2)+1
        Bloqade.plot!(ax1, Δ2_single_defect[idx])
    else
        Bloqade.plot!(ax2, Δ2_single_defect[idx])
    end
end

fig3

# We can confirm that the waveforms produce the desired domain walls for the LGT states, 
# by simulating the dynamics governed by the waveforms, followed by plotting their density profiles:

dens2 = get_average_rydberg_densities(Δ2_single_defect, Ω2)
dens3 = get_average_rydberg_densities(Δ2_two_defects, Ω2)

fig4 = Figure(size=(800, 320))
ax1 = Axis(fig4[1, 1], xlabel="site index", ylabel="Average Rydberg densities")
ax2 = Axis(fig4[2, 1], xlabel="site index", ylabel="Average Rydberg densities")
barplot!(ax1, 1:N, dens2[end])
barplot!(ax2, 1:N, dens3[end])
fig4

# Again, we see that the Rydberg density at the defects are not exactly zero, but the prepared states,
# as we shall see below, serve as good initial states to study the propagation of particle-antiparticle pairs in LGT. 

# We define the very last piece in the Rabi frequency and detuning that govern the time evolution of the Rydberg chain with defects:

Δq2 = -π ;
tq2 = 40/Ωmax ; 

Δ3_single_defect = map(1:length(atoms)) do idx
    append(Δ2_single_defect[idx], constant(duration=tq2, value=Δq2))
end       
Δ3_two_defects = map(1:length(atoms)) do idx
    append(Δ2_two_defects[idx], constant(duration=tq2, value=Δq2))
end       

Ω3 = append(Ω2, constant(duration=tq2, value=Ωmax));

# Again, as an example, for the case with a single defect, we show the detuning for the central site, 
# which is the defect, and those for other sites separately below:

fig5 = Figure(size=(960, 320))
ax1 = Axis(fig5[1, 1], xlabel="Time", ylabel="Detuning", xgridvisible=true, ygridvisible=true)
ax2 = Axis(fig5[1, 2], xlabel="Time", ylabel="Detuning", xgridvisible=true, ygridvisible=true)

Bloqade.plot!(ax1, Δ3_single_defect[floor(Int, N/2)+1], label="Detuning for the central site")
Bloqade.plot!(ax2, Δ3_single_defect[1], label="Detunings for other sites")
Bloqade.plot!(ax1, Ω3; label="Rabi frequency for all sites")
axislegend(ax1, position=:rb)
axislegend(ax2, position=:rb)
fig5

# ### Simulating Particle-Antiparticle Pairs in LGT Dynamics 

# With the waveforms defined above, we can run the simulation to evolve the Rydberg chains with defects and collect the final Rydberg densities:
densities_single_defect = get_average_rydberg_densities(Δ3_single_defect, Ω3);
densities_two_defects = get_average_rydberg_densities(Δ3_two_defects, Ω3);

D_single_defect = hcat(densities_single_defect...);
D_two_defects = hcat(densities_two_defects...);

# To better visualize the propagation of particle-antiparticle pairs, 
# we shall only show the Rydberg densities starting from the time point when the ground state of the defect chain is prepared:

ind0 = 3550;
D_single_defect = D_single_defect[:, ind0:end];
D_two_defects = D_two_defects[:, ind0:end];

clocks = 0:1e-3:Ω3.duration;
clocks = clocks[ind0: end];

# Then we plot the Rydberg density as a function of time, where the two panels correspond to the cases with single and two defects respectively:
yticks = range(clocks[1], stop=clocks[end], length=10);
yticks = [string(ytick)[1:4] for ytick in yticks][end:-1:1];

fig6 = Figure(size=(960, 640))
ax1 = Axis(fig6[1, 1], xlabel = "sites", ylabel = "time (μs)", yticks = (range(1, stop = length(clocks), length = 10), string.(round.(range(clocks[1], stop=clocks[end], length=10); digits=2))))
ax2 = Axis(fig6[1, 2], xlabel = "sites", yticklabelsvisible=false)
clims = (0.0, 1.0)
img1 = image!(ax1, D_single_defect[end:-1:1, :], colormap = :inferno, colorrange=clims, interpolate=false)
img2 = image!(ax2, D_two_defects[end:-1:1, :], colormap = :inferno, colorrange=clims, interpolate=false)
Colorbar(fig6[1, 3]; limits=clims)
fig6
# From the left panel, we can observe a light-cone-shaped region originating from the particle-antiparticle pair in the vacuum. 
# At the right panel, we show the interference of two light cones, which produces an additional change of periodicity corresponding to the elastic scattering. 
# When the particle or antiparticle reaches the boundary of the chain, it will be scattered back as observed. 
# For more details, the interested readers are referred to the paper [F. M. Surace et al. (10.1103/PhysRevX.10.021041)](https://journals.aps.org/prx/pdf/10.1103/PhysRevX.10.021041). 

# In summary, we have shown that the ground state and the dynamics of certain LGT can be simulated by a 1D chain of Rydberg atoms.
# More interestingly, defects can be introduced by locally addressing certain atoms in the chain, 
# and with that we can simulate the propagation of particle-antiparticle pairs in the LGT dynamics.


