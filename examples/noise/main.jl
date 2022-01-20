using EaRyd
using CairoMakie
using EaRydPlots

# In this example, we examine the adiabatic preparation of a spin Hamiltonian with noise parameters to show the erro propagation feature.  


# We first prepare the adiabatic pulse sequence as two piecewise linear functions
# the first part considers Hamiltonian parameters with time-dependent global noise (where each atom site is subject to the same noise)

# define the rabi waveform
Ω_src = 2.3 * 2 * pi
# the noise strength can be specified by \pm
Ω_max = Ω_src ± 0.001
Ω = piecewise_linear(clocks=[0.0, 0.252, 1.052, 1.6], values=[0.0, Ω_max , Ω_max , 0])

# define the noisy detuning waveform
U = Ω_src / 2.3 ± 0.001
Δ = piecewise_linear(clocks=[0.0, 0.252, 1.052, 1.6], values=[-6*U, -6*U, 2*U , 2*U])

# We prepare a square lattice of 9 atoms 
# create the atom positions
nx, ny = 3, 3
nsites = nx*ny

atoms = generate_sites(SquareLattice(), nx, ny, scale=9.629)


# We construct the Rydberg Hamiltonian from the defined noisy rabi and detuning waveforms
h = rydberg_h(atoms; C=2 * pi * 858386, Δ, Ω)

using SparseArrays
# We evolve the system from the zero state using the ODE solver to a final time t = 1.6 microseconds
prob = ODEEvolution(zero_state(9), 1.6±0.0, h)
emulate!(prob) # run the time evolution directly

# We compute the Rydberg probability for each site, where each of the results contains an error bar
densities = map(1:nsites) do i
    real(expect(put(nsites, i=>Op.n), prob.reg))
end



# the second part considers Hamiltonian parameters with time-dependent local noise (where each atom site is subject to independent noise)


# each atom is subject to an indenpendent noise
Ω = map(1:9) do idx
    piecewise_linear(clocks=[0.0, 0.252, 1.052, 1.6], values=[0.0, Ω_max, Ω_max , 0])
end

Δ = map(1:9) do idx
    piecewise_linear(clocks=[0.0, 0.252, 1.052, 1.6], values=[-6*U, -6*U, 2*U , 2*U])
end


# define the Hamiltonian 

h = rydberg_h(atoms; C=2 * pi * 858386, Δ, Ω)
prob = ODEEvolution(zero_state(9), 1.6±0.0, h)
emulate!(prob) # run the time evolution directly

# measure the Rydberg density for each site
densities = map(1:nsites) do i
    real(expect(put(nsites, i=>Op.n), prob.reg))
end
