# write your EaRyd example with Literate.jl here
using EaRyd
using CairoMakie
using EaRydPlots
using SparseArrays


# In this example, we examine the adiabatic preparation of a spin Hamiltonian with noise parameters to show the erro propagation feature.  


# We first prepare the adiabatic pulse sequence as two piecewise linear functions
# the first part considers Hamiltonian parameters with time-dependent global noise (where each atom site is subject to the same noise)

# define the rabi waveform
# the noise strength can be specified by \pm
total_time = 2.9
Ω_max = 2π * 4.3 ± 1
Ω = piecewise_linear(clocks=[0.0, 0.3, 2.6, total_time], values=[0.0, Ω_max , Ω_max , 0])


# define the noisy detuning waveform
U = 2π * 15.0 ± 8
Δ = piecewise_linear(clocks=[0.0, 0.3, 2.6, total_time], values=[-U, -U, U , U])

# We prepare a square lattice of 9 atoms 
# create the atom positions
nx, ny = 3, 3
nsites = nx*ny
ax= 6.7 
ay=6.7
dx= 0±0.1
dy=0±0.1


#atom_coordinate= [(dx,  dx)]

#for ii= 2: nx
 #   push!(atom_coordinate, (ax*(ii-1)+ dx, 0* ay+ dy))
#end

#for jj=2:ny
#for ii= 1: nx
 #   push!(atom_coordinate, (ax*(ii-1)+ dx, (jj-1)*ay+ dy))
#end
#end


#atoms = generate_sites(SquareLattice(), nx, ny, scale=6.7)
atoms= map(loc -> loc .± (0.1, 0.1), generate_sites(SquareLattice(), nx, ny, scale=6.7))



# We construct the Rydberg Hamiltonian from the defined noisy rabi and detuning waveforms
#h = rydberg_h(atom_coordinate; C=2 * pi * 858386, Δ, Ω)


h = rydberg_h(atoms; C=2 * pi * 858386, Δ, Ω)


reg = zero_state(9)
# We evolve the system from the zero state using the ODE solver to a final time t = 1.6 microseconds
prob = ODEEvolution(reg, total_time, h; dt=1e-3, adaptive=false)
emulate!(prob) # run the time evolution directly

# We compute the Rydberg probability for each site, where each of the results contains an error bar
densities = map(1:nsites) do i
    real(expect(put(nsites, i=>Op.n), prob.reg))
end


print(densities)
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
print(densities)