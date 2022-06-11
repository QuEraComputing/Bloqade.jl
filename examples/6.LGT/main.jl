# # Simulating lattice gauge theory with Rydberg atoms
# ## Introduction
# In the [Adiabatic Evolution](@ref adiabatic-evolution) page, we show how to prepare ordered
# ground states in the Rydberg system. In the [@ref quantum-scar] page, we see that the 
# evolution of an ordered ground state exhibits revival behavior, which oscillates between two 
# patterns of Rydberg densities. 
#
# In this tutorial, we shall simulate the 1D latticee gauge theory (LGT) with a 1D Rydberg atom chain. In particular, we shall illustrate that the ground state of the studied LGT corresponds to the ordered state of the Rydberg chain, and the quantum scar is nothing but the string-inversino mechanism in the gauge theories. More interestingly, via creating local defects in the chain, we can simulate the propagation of particle-antiparticle pairs. This tutorial is inspired by the paper https://journals.aps.org/prx/pdf/10.1103/PhysRevX.10.021041

# We first import the required packages 

using Bloqade
using PythonCall
# using PyCall
using StatsBase
using Distributed
using BitBasis

# plt = pyimport("matplotlib.pyplot");


using Dates
start = time()
using Plots


# # Prepare the ground state
# The prerequisite of this experiment is to prepare the so-called vaccum state for the LGT. In the language of Ising model, it is nothing but the anti-ferromagnetic state. We know how to prepare the 1D AF state from the 0th use case, and starting from the ground state, we simulate the dynamics of the LGT using the Rydberg system.

# More background of mapping LGT to Rydberg....

# ## Prepare the lattice
# In order to simulate the dynamics of the LGT on Rydberg system, we first prepare a 1D lattice with 21 sites

a = 5.5 # 7.0, 5.5
N = 21 # 13 # 17 # needs to be 4n+1

atoms = generate_sites(ChainLattice(), N, scale=a)
subspace = blockade_subspace(atoms, a)

### good paras
total_time = 3.5; 
Ωmax = 2π * 5;
Δ1 = -2π * 10;
Δ2 = 2π * 10;

dt = 1e-3 # μs

Ωq = Ωmax
tq = pi/Ωq

Ωq2 = Ωmax
tq2 = 40/Ωq2

Δ1 = piecewise_linear(clocks = [0.0, 0.2, total_time], values = [Δ1, Δ1, Δ2]); 
Ω1 = piecewise_linear(clocks = [0.0, 0.2, total_time-0.0001, total_time], values = [0.0, Ωmax, Ωmax, 0]);

# The profiles of the waveforms are shown below

# fig1, ax = plt.subplots(1, 1, figsize = (10,4))
# Bloqade.plot!(ax, Ω1)
# Bloqade.plot!(ax, Δ1)
# ax.grid()

# ax.legend(["Ω", "Δ"])

# fig1

# We can confirm that the waveforms produce the desired initial state of the LGT, by simulating the dynamics governed by the LGT Hamiltonian, followed by plotting the density profile. 

function Get_Ryd_density(Δ, Ω)
    h = rydberg_h(atoms; Δ=Δ, Ω=Ω)
    reg = zero_state(subspace)

    # We can then simulate the time evolution of the quantum state using an ODE solver:
    duration = Δ.duration
    # prob = SchrodingerProblem(reg, duration, h, adaptive=false, dt=1e-5, progress=true); # if adaptive bigger than dt, then error 
    prob = SchrodingerProblem(reg, duration, h, progress=true); # if adaptive bigger than dt, then error 
    integrator = init(prob, Vern8());
    # Then, we measure the real-time expectation value of the Rydberg density and entanglement entropy: 

    densities = []
    for _ in TimeChoiceIterator(integrator, 0.0:dt:duration)
        normalize!(reg) # can normalize 
        push!(densities, rydberg_density(reg)) 
    end

    # print(norm(reg))

    # bitstring_hist(reg; nlargest = 20)
    return reg, densities

    # check reg norm, norm(reg)

end

reg1, dens1 = Get_Ryd_density(Δ1, Ω1)
# bitstring_hist(reg1; nlargest = 20)
# fig1, ax = plt.subplots(figsize = (10, 4))
# ax.bar(1:N, dens1[end])
# fig1


# # Preparation of ground state of Rydberg chain with defect

# Up to now, what we have done is simply reproducing the example in here. 

# Following the above section, next we prepare the ground state with a defect at the center. To achieve that, we need the second part of the pulse program, where we will need the local control for the center atom. More specifically, we shall apply a -pulse for the central atom such that it transition from the Rydberg state to ground state.

# ## Define the site-dependent waveform

Δ2 = append(Δ1, constant(duration=tq, value=0))
Ω2_single_defect = map(1:length(atoms)) do idx
    if idx == floor(Int, N/2)+1
        append(Ω1, constant(duration=tq, value=Ωq))
    else
        append(Ω1, constant(duration=tq, value=0))
    end
end

Ω2_two_defects = map(1:length(atoms)) do idx
    if idx == floor(Int, N/3) || idx == floor(Int, N-N/3)+1
        append(Ω1, constant(duration=tq, value=Ωq))
    else
        append(Ω1, constant(duration=tq, value=0))
    end
end


# The profile of the site-dependent Rabi frequencies are shown below

# fig2, ax = plt.subplots(1, 1, figsize = (10,4))
# for Ωp in Ω2_single_defect
#     Bloqade.plot!(ax, Ωp)
# end

# ax.grid()

# ax.legend(range(1,length(Ω2_single_defect)))

# fig2

# We can confirm that the waveforms produce the desired initial state of the LGT, by simulating the dynamics governed by the LGT Hamiltonian, followed by plotting the density profile. 

reg2 = Get_Ryd_density(Δ2, Ω2_single_defect)[1]
bitstring_hist(reg2; nlargest = 20)

# # Propagation of particle-antiparticle pairs
# Some physics background
#
# Here, we create a defect in the Rydberg array, and prepare the system into the ground state, and observe the propagation of the particle and anti-particle pairs.


# ## Define the site-dependent waveform


Δ3 = append(Δ2, constant(duration=tq2, value=-π))
Ω3_single_defect = map(1:length(atoms)) do idx
    append(Ω2_single_defect[idx], constant(duration=tq2, value=Ωq2))
end       
Ω3_two_defects = map(1:length(atoms)) do idx
    append(Ω2_two_defects[idx], constant(duration=tq2, value=Ωq2))
end       


# The profile of the site-dependent Rabi frequencies are shown below


# fig3, ax = plt.subplots(1, 1, figsize = (10,4))
# for Ωp in Ω3_single_defect
#     Bloqade.plot!(ax, Ωp)
# end

# ax.grid()

# ax.legend(range(1,length(Ω3_single_defect)))

# fig3

# ## Simulation particle-antiparticle pairs in LGT dynamics 

# # h = rydberg_h(atoms; Δ = Δ3, Ω = Ω3)
# h = rydberg_h(atoms; Δ = Δ2, Ω = Ω2)
# reg = zero_state(subspace)
# prob = SchrodingerProblem(reg, Δ3.duration, h);
# integrator = init(prob, Vern8());

# densities = []
# for _ in TimeChoiceIterator(integrator, 0.0:dt:Δ3.duration)
#     push!(densities, rydberg_density(reg))
# end

densities_single_defect = Get_Ryd_density(Δ3, Ω3_single_defect)[2]
densities_two_defects = Get_Ryd_density(Δ3, Ω3_two_defects)[2]

# To better visualize the propagation of particle-antiparticle pairs, we set the initial time to be when the ground state of the defect chain is prepared. 

final = time()
print(final-start)

ind0 = 3550 # 3500

D_single_defect = hcat(densities_single_defect...)
D_two_defects = hcat(densities_two_defects...)

clocks = 0:1e-3:Δ3.duration

D_single_defect = D_single_defect[:, ind0:end]
D_two_defects = D_two_defects[:, ind0:end]
clocks = clocks[ind0: end]
# D = D[:, 1: ind0]
# clocks = clocks[1: ind0]

# Then we plot the Rydberg density of the time as a function of time

hms = [heatmap(real(transpose(D_single_defect)), legend=:none), heatmap(real(transpose(D_two_defects)), legend=:none)] # Make heatmaps without legends
l = @layout[grid(1,2) a{0.05w}] # Stack a layout that rightmost one is for color bar

p = plot(hms..., heatmap((0:0.01:1).*ones(101,1), legend=:none, xticks=:none, yticks=(1:10:101, string.(0:0.1:1))), layout=l) # Plot them set y values of color bar accordingly

# julia> savefig(p, "plot.png")

# fig4, ax = plt.subplots(1, 2, figsize = (6, 3))

# ax[0].imshow(real(transpose(D_single_defect)[end:-1:1, :]), interpolation = "nearest", aspect = "auto", extent = [0.5, N + 0.5, 0, Δ3.duration])
# ax[0].set_ylabel("time (μs)")
# ax[0].set_xlabel("site")
# # ax.set_yticks(0:0.4:Δ3.duration)
# ax[0].set_xticks(1:N)

# shw = ax[1].imshow(real(transpose(D_two_defects)[end:-1:1, :]), interpolation = "nearest", aspect = "auto", extent = [0.5, N + 0.5, 0, Δ3.duration])
# # ax[1].set_ylabel("time (μs)")
# ax[1].set_yticks([])
# ax[1].set_xlabel("site")
# # ax.set_yticks(0:0.4:Δ3.duration)
# ax[1].set_xticks(1:N)


# # bar = fig4.colorbar(shw)
# fig4

# As we see, the Rydberg chain can simulate the particle-antiparticle pairs in LGT dynamics by locally detuning one or a few atoms to simulate the defect. 

