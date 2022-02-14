# # Background

using EaRyd
using CairoMakie
using EaRydPlots


# # Two dimensional case 


# In this example, we examine the adiabatic preparation of a spin Hamiltonian.  

# We first prepare the adiabatic pulse sequence as two piecewise linear functions
# define the rabi waveform
total_time = 2.9
Ω_max = 2π * 4.3
Ω = piecewise_linear(clocks=[0.0, 0.3, 2.6, total_time], values=[0.0, Ω_max , Ω_max , 0])
draw(Ω)

# define the detuning waveform
U = 2π * 15.0
Δ = piecewise_linear(clocks=[0.0, 0.3, 2.6, total_time], values=[-U, -U, U , U])
draw(Δ)


# We prepare a square lattice of 9 atoms 
# create the atom positions
nx, ny = 3, 3
nsites = nx*ny
atoms = generate_sites(SquareLattice(), nx, ny, scale=6.7)


# We construct the Rydberg Hamiltonian from the defined rabi and detuning waveforms
h = rydberg_h(atoms; Δ, Ω)

# we now declare our register, a register is an object stores the information about the initial quantum state
reg = zero_state(9)


# We evolve the system from the zero state using the ODE solver to a final time t = 1.6 microseconds
prob = ODEEvolution(reg, total_time, h; dt=1e-3, adaptive=false)

# now we can run the evolution, if you don't need to measure anything during the evolution, you can just run
# `emulate!(prob)`
#but here we will be plotting the total number of rydberg states at each time

densities = []
for info in prob
    push!(densities, [expect(put(nsites, i=>Op.n), info.reg) for i in 1:nsites])
end
D = hcat(densities...)

# now we can plot the Rydberg density as a function of time

clocks = [t for t in 0:1e-3:total_time]
heatmap(clocks, 1:9, D'; axis=(xlabel="iterations", ylabel="rydberg density per site"))

# plot the histogram of the most frequent results from sampling the final state

bitstring_histgram(prob.reg; nlargest=20)




# # One dimensional case, 
# maybe this comes first before the 2D case

# this is the Rabi pulse
total_time = 5
Ω_max = 2π * 2
Ω = piecewise_linear(clocks=[0.0, 0.1, 3.1, 3.2, total_time], values=[0.0, Ω_max , Ω_max , 0, 0])
draw(Ω)

# this is the Detuning pulse 
U1 = -2π * 4
U2 = 2π * 9
Δ = piecewise_linear(clocks=[0.0, 0.6, 3.2, total_time], values=[U1, U1, U2 , U2])
draw(Δ)

# We prepare a chain lattice of 11 atoms 
# create the atom positions, scale = 5.72 for Z2, 3.57 for Z3, 2.87 for Z4
nsites = 9
atoms = generate_sites(ChainLattice(), nsites, scale=5.72)
h = rydberg_h(atoms; Δ, Ω)
reg = zero_state(9)
prob = ODEEvolution(reg, total_time, h; dt=1e-3, adaptive=false)


densities = []
for info in prob
    push!(densities, [expect(put(nsites, i=>Op.n), info.reg) for i in 1:nsites])
end
D = hcat(densities...)


clocks = [t for t in 0:1e-3:total_time]
heatmap(clocks, 1:9, D'; axis=(xlabel="iterations", ylabel="rydberg density per site"))


# we can also plot the bitstring
bitstring_histgram(prob.reg; nlargest=20)




