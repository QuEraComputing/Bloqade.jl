using Bloqade

# Simulate nonequilibrium dynamics of a 12-site ring of Rydberg atoms.
#  This script generates Figure 2 in the example RydbergBlockade.

# First, define the geometry of the system, a ring of nsites.
nsites = 14;    # number of sites in the chain
distance = 9    # Distance between atoms, in microns

R = distance / (2 * sin(2 * pi / (nsites) / 2))                                       # Radius of the circle, using a little trigonometry
pos = [(R * sin(i * 2 * pi / (nsites)), R * cos(i * 2 * pi / (nsites))) for i in 1:nsites] # Positions of each atom
atoms = BloqadeLattices.AtomList(pos)                                         # Define the atom positions as an AtomList.

# Define the pulse by specifying a maximum omega and delta.
# Note that these parameters are computed using pulser_compare.py
Omega_max = 9.179086064116243
delta_0 = -20.397969031369428
delta_f = 9.179086064116243

# Define the timescales of the pulse
Trise = 2.0
Tramp = 6.0
Tfall = 2.0
total_time = Trise + Tramp + Tfall

# Define the detuning and Rabi fields as a function of time
Δ = piecewise_linear(clocks = [0.0, Trise, Trise + Tramp, total_time], values = [delta_0, delta_0, delta_f, delta_f]);
Ω = piecewise_linear(clocks = [0.0, Trise, Trise + Tramp, total_time], values = [0, Omega_max, Omega_max, 0]);

# Define the Hamiltonian and problem instance
h = rydberg_h(atoms; C = 5420158.53, Δ, Ω)
reg = zero_state(nsites);
problem = SchrodingerProblem(reg, total_time, h);

# Evolve the state to the final state
@time begin
    emulate!(problem)
end

# Compute the wavefunction probabilities to cross-check with pulser and quspin code
probability = abs.(reg.state) .^ 2
for p in reverse(sort(probability, dims = 1))[1:10]
    print(p, '\n')
end
