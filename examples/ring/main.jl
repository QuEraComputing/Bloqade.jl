# # Ring Emulation
using Bloqade

# Use matplotlib to generate plots
# for windows users you may need to use the following commands
# ```julia
# matplotlib = pyimport("matplotlib")
# matplotlib.use("TkAgg")
# ```

using PythonCall
plt = pyimport("matplotlib.pyplot")
plt.rcParams["font.size"] = 22

# Simulate nonequilibrium dynamics of a 12-site ring of Rydberg atoms.
# This script generates Figure 2 in the example RydbergBlockade.


# First, define the geometry of the system, a ring of 12 sites.
nsites = 12;    # 12 site chain
distance = 7    # Distance between atoms, in microns

R = distance/(2*sin(2*pi/(nsites)/2))                                       # Radius of the circle, using a little trigonometry
pos = [(R*sin(i*2*pi/(nsites)), R*cos(i*2*pi/(nsites)) ) for i in 1:nsites] # Positions of each atom
atoms = AtomList(pos)                                                       # Define the atom positions as an AtomList.

# Build the Hamiltonian by defining the atom list and the Rabi frequency
# Here, the Rabi frequency Ω is π, such that one oscillation occurs in 0.5usec.
h = rydberg_h(atoms;C = 2π * 858386, Ω=π)


# # Emulate the problem

# State in the full space
init_state = zero_state(nsites)
# State in the blockade subspace
space = blockade_subspace(atoms,distance*1.1)   # Compute the blockade subspace
init_state2 = zero_state(space)                       # Define the initial state in the blockade subspace.


# Define the time steps
Tmax = 10.
nsteps = 5001
times = LinRange(0,Tmax,nsteps)


# Time evolve the system in the full space
prob = SchrodingerProblem(init_state, Tmax, h, dt = Tmax/(nsteps-1) , adaptive = false);
integrator = init(prob, Vern6());

densities = []
for _ in TimeChoiceIterator(integrator, 0.0:Tmax/(nsteps-1):Tmax)
    push!(densities, rydberg_density(init_state, 1))
end


# Time evolve the system in the subspace
prob2 = SchrodingerProblem(init_state2, Tmax, h, dt = Tmax/(nsteps-1) , adaptive = false);
integrator2 = init(prob2, Vern8());

densities2 = []
for _ in TimeChoiceIterator(integrator2, 0.0:Tmax/(nsteps-1):Tmax)
    push!(densities2, rydberg_density(init_state2, 1))
end

# Plot the data
fig, ax = plt.subplots(figsize=(8,6))

ax.plot(times,real(densities),"k",label="Full space")
ax.plot(times,real(densities2),"r--",label="Subspace")
ax.axis([0,Tmax,0,0.45])
ax.set_xlabel("Time (us)")
ax.set_ylabel("Rydberg density")
plt.tight_layout()
ax.legend()

inset_axes = pyimport("mpl_toolkits.axes_grid1.inset_locator")
ax2 = inset_axes.inset_axes(ax,width="20%",height="30%",loc="lower right",borderpad=1)
ax2.plot(times,real(densities - densities2))
ax2.axis([0,0.6,-0.0008,0.0008])
ax2.set_ylabel("Difference",fontsize=12)
ax2.set_yticks(LinRange(-0.0008,0.0008,5),fontsize=12)
ax2.set_xticks([0,0.2,0.4,0.6],fontsize=12)
fig
