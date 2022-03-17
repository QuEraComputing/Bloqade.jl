using EaRyd

# Simulate nonequilibrium dynamics of a 12-site ring of Rydberg atoms.
#  This script generates Figure 2 in the example RydbergBlockade.


# First, define the geometry of the system, a ring of 12 sites.
nsites = 12;    # 12 site chain
distance = 7    # Distance between atoms, in microns

R = distance/(2*sin(2*pi/(nsites)/2))                                       # Radius of the circle, using a little trigonometry
pos = [(R*sin(i*2*pi/(nsites)), R*cos(i*2*pi/(nsites)) ) for i in 1:nsites] # Positions of each atom
atoms = EaRydLattices.AtomList(pos)                                         # Define the atom positions as an AtomList.

# Build the Hamiltonian by defining the atom list and the Rabi frequency
# Here, the Rabi frequency Ω is π, such that one oscillation occurs in 0.5usec.
h = rydberg_h(atoms;C = 2π * 858386, Ω=π)


# # Emulate the problem

# Define the initial state in the full basis
init = zero_state(nsites)

# Do the same thing in the blockade subspace...
space = blockade_subspace(atoms,distance*1.1)   # Compute the blockade subspace
init2 = zero_state(space)                       # Define the initial state in the blockade subspace.


# 
Tmax = 10
nsteps = 20
iteration = 1:nsteps
ts = [Tmax/nsteps for _ in iteration];
hs = [h for _ in iteration];
clocks = cumsum(ts);


prob = KrylovEvolution(init, ts, hs)
prob2 = KrylovEvolution(init2, ts, hs)


# Then we measure the real-time expectation value of Rydberg density, domain wall density, and entanglement entropy. 
# These data are stored in the matrix or vector below. 

data_out = zeros(3, length(iteration))
data_out[1,:] = clocks
@time begin
    for info in prob
        data_out[2, info.step] = expect(put(nsites, 1=>Op.n), info.reg)
    end
end

@time begin
    for info in prob2
        data_out[3, info.step] = expect(put(nsites, 1=>Op.n), info.reg)
    end
end


using DelimitedFiles

fil = open("density_matrix_data.txt","w")
writedlm(fil,data_out)
close(fil)
