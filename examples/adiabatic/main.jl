using EaRyd
using CairoMakie
using EaRydPlots

# In this example, we examine the adiabatic preparation of a spin Hamiltonian.  


# We first prepare the adiabatic pulse sequence as two piecewise linear functions
# define the rabi waveform
Ω_max = 2.3 * 2 * pi
Ω = piecewise_linear(clocks=[0.0, 0.252, 1.052, 1.6], values=[0.0, Ω_max , Ω_max , 0])

# define the detuning waveform
U = Ω_max / 2.3
Δ = piecewise_linear(clocks=[0.0, 0.252, 1.052, 1.6], values=[-6*U, -6*U, 2*U , 2*U])


# We prepare a square lattice of 9 atoms 
# create the atom positions
nx, ny = 3, 3
nsites = nx*ny

atoms = generate_sites(SquareLattice(), nx, ny, scale=9.629)


# We construct the Rydberg Hamiltonian from the defined rabi and detuning waveforms
h = rydberg_h(atoms; C=2 * pi * 858386, Δ, Ω)

# We evolve the system from the zero state using the ODE solver to a final time t = 1.6 microseconds
prob = ODEEvolution(zero_state(9), 1.6, h)
emulate!(prob) # run the time evolution directly

# We compute the Rydberg probability and plot the histogram of the most frequent results from sampling the final state 
densities = map(1:nsites) do i
    real(expect(put(nsites, i=>Op.n), prob.reg))
end

bitstring_histgram(prob.reg; nlargest=20)
