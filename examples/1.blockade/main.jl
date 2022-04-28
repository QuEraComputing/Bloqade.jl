# # [Rydberg Blockade](@id blockade)

# The Rydberg blockade mechanism is a crucial component of the operation of neutral atom computers. By including a term in the effective Hamiltonian which adds an interaction between adjacent Rydberg atoms, we may construct gates and nontrivial dynamics which create entanglement and correlation across the system.

# The term is a van der Waals interaction where there is an energy shift if two adjacent atoms are in the Rydberg state. The Van der Waals interaction is

# ```math
# V_{ij} = \frac{C_6}{|\vec r_i - \vec r_j|^6}\hat n_i \hat n_j
# ```

# where ``\hat n_i=|r_i\rangle\langle r_i|`` is the number operator on the ``i``th site, which is 1 if the atom is in the Rydberg state, and zero if the atom is in the ground state. The coefficient ``C_6 = 2\pi\times 862690 MHz/\mu m^6`` is the interaction strength; characteristically, this interaction has a strength ``\approx 10MHz`` for two atoms separated by ``10\mu m``, a similar scale to the Rabi coupling between the ground and Rydberg state. Crucially, this can be seen as an energy shift on atom ``j``, conditional on the state of atom ``i``, and so can be used, in a loose sense, as a conditional logical gate. This is because the adjacent atom forces the atom to be off-resonant with the laser field.

# This conditional drive can be seen given the following dynamics. Suppose two atoms are close to each other (``< 10 \mu m``) and so interact under van der Waals. The left atom is either in a Rydberg state, or in the ground state, and the right atom is originally in the ground state. Then, a Rabi drive is applied to the right atom, which couples the atom's ground state to the Rydberg state. For this example, we choose a Rabi drive of ``\Omega=2\pi\times 0.5``MHz and distance between atoms ``|\vec r_i- \vec r_j| = 7\mu m``, which gives a conditional detuning of ``\approx 50``MHz. When the left atom is in the ground state (black, top), there are no interactions and the state of the right atom oscillates between the ground state and Rydberg state; for a particular choice of timing this executes a ``\pi`` pulse, flipping the right atom from the ground to Rydberg state. However, when the left atom is in the Rydberg state (red, bottom), there is a large detuning on the right atom, which causes the transfer to the Rydberg state to be strongly suppressed. In this case, the right atom (up to perturbative corrections) is never in the Rydberg state.



# ![RydbergBlockade](../../../assets/RydbergBlockade.png)

# This conditional energy shift is the basis of the **Rydberg Blockade**. Because of the large energy shift from having two adjacent atoms in the Rydberg state, evolution from an atomic ground state with local Rabi couplings between ground and Rydberg is restricted to a low energy subspace of states where no two adjacent atoms are in the Rydberg state. Furthermore, because the interaction strength with distance is a large power law, one can define a characteristic scale set by the Rabi coupling. If two atoms are close such that the conditional detuning is much larger than the Rabi coupling, one can consider the atoms to be blockading each other, and both atoms cannot simultaneously be in the Rydberg state. In contrast, if two atoms are far away, the two atoms never blockade each other and both atoms can simultaneously be in the Rydberg state.

# ## Energy truncation and Blockade subspace

# ![EnergyTruncation](../../../assets/bloqade_subspace.png)

# In this way, the low energy classical states of the Rydberg Hamiltonian (``\Omega=0``, ``\Delta=0``) for a given array of atoms are **independent sets** of a **unit disk graph** defined by the positions of the atoms. A unit disk graph is a set of vertices and edges, where vertices represent every atom, and there are edges if the distance between vertices is less than some unit disk radius ``|\vec r_i- \vec r_j|<R_u``. The lowest energy states are representative of independent sets, where Rydberg excitations are in the independent set and no two Rydberg excitations are within some radius. The second energy band are sets with a single independent set violation, where there is equivalently two Rydberg excitations within the unit disk radius of each other. Higher and higher bands represent more and more independent set violations. Note that this band structure is dependent on the arrangement of atoms, and for arbitrary configurations this band structure may not be clear.

# Instead of doing quantum dynamics in the full ``2^N`` Hilbert space, we may take advantage of the energy structure of the classical Hamiltonian to reduce the computational difficulty. The simplest scheme is to truncate the Hilbert space to the low energy subspace, and exclude all states above a certain energy. Given the natural band structure of the classical Hamiltonian, we may simply truncate the Hilbert space to the subspace of independent sets of the unit disk graph. Equivalently, this is the **blockade subspace**, where atoms within the blockade radius are excluded from both being in the Rydberg state.




# The validity of the energy truncation subspace is governed by the strength of off-diagonal matrix elements coupling the low energy subspace to the high energy one. For the Rydberg Hamiltonian, these off-diagonal elements ``|1r\rangle\leftrightarrow|rr\rangle`` have a strength ``\Omega``. In order to preserve dynamics within the subspace, the energy difference between states within the blockade subspace (e.g., ``|1r\rangle``) and outside (``|rr\rangle``) must be much larger than the Rabi strength. Formally,

# ```math
# \Omega \ll \frac{C_6}{R_\text{min}^6}
# ```
# where ``R_\text{min}`` is the minimum the unit disk radius can be without removing vertices to the graph, or equivalently the maximum distance between any two vertices that are connected. As long as this condition holds, the exact dynamics in the full Hilbert space should be closely approximated by the approximate dynamics in the blockade subspace, as the mixing terms only couple to low energy states.

# ## ``\Omega_b``, unit disk radius, and blockade radius

# While in principle the unit disk radius can be made equal to the blockade radius, the two values should be differentiated. For any two atoms within the unit disk radius ``R_u``, the energy scale of having both in the Rydberg state must be much larger than ``\Omega``, as described above. However, for any two atoms outside of the unit disk radius, the energy scale of having both be in the Rydberg state must be much smaller than ``\Omega``, as the two atoms should not be blockaded. This condition sets a lower limit on the value of ``\Omega``, as there are still ``1/R^6`` interactions between nearby atoms which may "accidentally" blockade each other if ``\Omega`` is too small. This lower limit guarantees that dynamics occur within the correct independent set subspace and is not affected by long range "Rydberg tails" which cause each independent set state to have a slightly different energy. The lower limit is set by ``R_\text{max}``, which is the maximum unit disk can be without adding any edges. This is equivalently the minimum distance between atoms that are _not_ within the unit disk radius:

# ```math
# \Omega \gg \frac{C_6}{R_\text{max}^6}
# ```

# Crucially, this lower limit is not necessary if one only cares about the validity of the Hilbert space truncation to the low energy space. This condition is additional, requiring that the dynamics within the low energy subspace effectively are in a degenerate independent set basis, and is not affected by long range ``1/R^6`` "Rydberg tails". One can than choose a value of ``\Omega`` which minimizes the error from both sources. Setting the ratio of perturbative errors equal to one, the optimal value of ``\Omega`` is set by

# ```math
# \Omega_b \equiv \frac{C_6}{R_b^6}\qquad\qquad\text{where}\qquad R_b = \sqrt{R_\text{min}R_\text{max}}
# ```

# Given that ``\Omega`` is usually fixed by the experemental system, it is usually more natural to set the unit disk radius by the Rabi drive

# ```math
# R_u = \bigg(\frac{C_6}{\Omega}\bigg)^{1/6} \sqrt{\frac{R_\text{min}}{R_\text{max}}}.
# ```

# There are now several characteristic distance scales:

# ``R_\text{min}`` is the minimum the unit disk radius can be without removing any edges from the unit disk graph.

# ``R_u`` is the unit disk radius of the graph. Usually, this is taken to be ``R_u=R_\text{min}``, though generally ``R_\text{min}\leq R_u<R_\text{max}``.

# ``R_b`` is the **blockade radius** of the system. This is the geometric mean ``R_b = \sqrt{R_\text{min}R_\text{max}}``

# ``R_\text{max}`` is the maximum the unit disk radius can be without adding any edges from the unit disk graph.

# ![BlockadRadius](../../../assets/bloqade_subspace_UDGradius.png)

# The ratio of the unit disk radius (dark red) to the blockade radius (red dashed) for several choices of atoms. Preferably, ``R_\text{max}/R_u\gg 1``, as large values mean that both perturbative conditions are preserved. For a 1d nearest neighbor line, if we choose  ``R_u=R_\text{min}``to be the lattice constant ``a`` and ``R_{max} =2a``, then the ratio between them is ``2``, and so ``(C_6/(2R_u)^6)\ll\Omega_b = (C_6/R_u^6)/8\ll (C_6/R_u^6)``. However, for an example arbitrary graph and if we choose the ``R_\text{min}`` and ``R_{max}`` to be the nearest- and next-nearest-neighbour atom distances respectively, this ratio is only ``1.15``, and so the perturbative limit ``0.42740\ll 0.6538\ll 1`` is not well-preserved. For this reason, there must be some care for choosing which graphs have appropriate dynamics within an approximately degenerate independent set subspace. Some graphs, like a 1d chain, are deep within the perturbative limits and so it is reasonable to expect that dynamics can ignore Rydberg tails. However, arbitrary graphs which have vertices close to the unit disk threshold may be sensitive to ``1/R^6`` Rydberg tails and may not have the expected dynamics within a degenerate independent set subspace.


# It should be noted and emphasized once again that this lower limit is only necessary if one needs the extra approximation of ignoring Rydberg interactions outside of the unit disk radius. If the subspace is chosen simply as an efficient way to speed up computation by truncating in energy, the value of ``\Omega`` does not have a lower bound. For example, if one wished to do dynamics of a ``\mathbb{Z}_3`` state, e.g. a next-nearest-neighbor chain, one can choose the blockade radius of a single site to truncate very high energy states, while choosing a value of ``\Omega`` which matches the blockade radius of the NNN chain.

# On the other hand,  if the user knows the value of ``\Omega``, then the corresponding blockade radius ``R_b`` can be computed by definition. Then the user can choose the unit disk radius ``R_u`` which is smaller than ``R_b`` and  satisfies the relation ``\Omega \ll \frac{C_6}{R_u^6}``. Then the user is able to compute the value of ``R_{max}`` by using the its relation with ``R_b`` and ``R_u``. This ensures the relation ``\Omega \gg \frac{C_6}{R_\text{max}^6}`` is satisfied.  Finally, the user is recommended to choose the subspace radius ``R_s`` to be ``R_u``, which makes the Hilbert space truncation only occurs for Rydberg excitations within the smallest radius. Thus accurate emulations are expected. For more information about setting the subspace radius, please refer to the section [subspace](@ref). 

# ## Example dynamics in the blockade subspace

# To emphisize the effectiveness of this independent set subspace, some example nonequilibrium dynamics are shown below, for a ring of 12 atoms seperated by ``7\mu m``. The minimum distance of atoms not within the blockade radius is ``\approx13.5230\mu m``, so that the the blockade radius is ``R_b=6.9\mu m\times \sqrt{13.32/6.9}\approx 9.58``. This sets the blockade energy scale is set to be ``\Omega_b \approx 1.01``MHz, and the perturbative limits to be ``\Omega_b\ll 7.29``MHz and ``\Omega_b\gg 0.14``MHz. This set of atoms can be defined by

using Bloqade
nsites = 12;    # 12 site chain
unit_disk_radius = 6.9    # Distance between atoms, in microns

R = unit_disk_radius/(2*sin(2*pi/(nsites)/2))                               # Radius of the circle, using a little trigonometry
pos = [(R*sin(i*2*pi/(nsites)), R*cos(i*2*pi/(nsites)) ) for i in 1:nsites] # Positions of each atom
atoms = AtomList(pos);                                                      # Define the atom positions as an AtomList.

blockade_radius = unit_disk_radius*sqrt(2);  # The blockade radius is the unit disk radius scaled by ~root2, the geometric factor for a line. Note that in principle this is slightly less due to the line actually being a finite-radius circle.
 
# The system is driven by a constant Rabi drive, which couples each atom's ground and Rydberg state. The value of ``\Omega\approx 1``MHz is set by the distance between vertices and is well within the perturbative limits set by the subspace, as ``1``MHz``\ll 7.296``MHz and so the blockade subspace should approximate exact dynamics faithfully. The Hamiltonian can be defined in bloqade with

Omega = 862690 / blockade_radius^6
h = rydberg_h(atoms;C = 2π * 862690, Ω=2*π*Omega);


# The system is initialized into the ground state of all atoms, which is the lowest energy of the classical Hamiltonian. We have two choices of basis: the first choice is the full Hilbert space of ``2^{12}`` elements, wheras the second basis is the blockade subspace, which excludes Rydberg excitations within the **subspace radius**. In principle, the subspace radius can be taken to be any value less than the blockade radius. For a subspace radius of zero, no states are excluded and one recovers exact dynamics. For a subspace radius anywhere between the blockade radius and ``R_\text{max}``, the subspace is the same, as there are no vertices within those radii. One reasonable value of the subspace radius for the ring of atoms is to simply choose it to be the unit disk radius. However, for more general graphs it may be reasonable to choose the subspace radius to be smaller than the unit disk radius and include extra states to improve the fidelity of the energy truncation. For example, for the next nearest neighbor line, it may be reasonable to choose the subspace radius to be half the unit disk radius, which includes high-energy NNN blockade states to improve numerical accuracy. See the [subspace](@ref) page for more details.
 
subspace_radius = unit_disk_radius

init_state = zero_state(nsites)                 # Define the initial state in the full space.
space = blockade_subspace(atoms,subspace_radius)# Compute the blockade subspace.
init_state2 = zero_state(space);                # Define the initial state in the blockade subspace.


# The blockade subspace has ``D=322`` elements, which means that computation is much faster. If the atoms were far apart and non-interacting, each atom would oscillate completely between its ground state and Rydberg state with a period of ``0.5 \mu``s. However, because adjacent atoms shift to the Rydberg state concurrently, they are dynamically blockaded, causing the maximum Rydberg density to only be 1/2, corresponding to an antiferromagnetic ``Z_2`` state. Note that because the ring has a translation symmetry, the Rydberg density is equal on all sites.


Tmax = 6.
nsteps = 2001
times = LinRange(0,Tmax,nsteps)



prob = SchrodingerProblem(init_state, Tmax, h, dt = Tmax/(nsteps-1) , adaptive = false);
integrator = init(prob, Vern6());

densities = [] # Time evolve the system in the full space
for _ in TimeChoiceIterator(integrator, 0.0:Tmax/(nsteps-1):Tmax)
    push!(densities, rydberg_density(init_state, 1))
end



prob2 = SchrodingerProblem(init_state2, Tmax, h, dt = Tmax/(nsteps-1) , adaptive = false);
integrator2 = init(prob2, Vern8());

densities2 = [] # Time evolve the system in the subspace
for _ in TimeChoiceIterator(integrator2, 0.0:Tmax/(nsteps-1):Tmax)
    push!(densities2, rydberg_density(init_state2, 1))
end

using PythonCall # Use matplotlib to generate plots
matplotlib = pyimport("matplotlib")
#matplotlib.use("TkAgg")
plt = pyimport("matplotlib.pyplot")

#Plot the data
#fig, ax = plt.figure(figsize=(8,6))

ax  = plt.subplot(1,1,1)
plt.plot(times,real(densities),"k",label="Full space")
plt.plot(times,real(densities2),"r--",label="Subspace")
ax.axis([0,Tmax,0,0.45])
plt.xlabel("Time (us)")
plt.ylabel("Rydberg density")
plt.tight_layout()
plt.legend()

inset_axes = pyimport("mpl_toolkits.axes_grid1.inset_locator")
ax2 = inset_axes.inset_axes(ax,width="20%",height="30%",loc="lower right",borderpad=1)
plt.plot(times,real(densities - densities2))
plt.axis([0,0.5,-0.001,0.003])
plt.ylabel("Difference",fontsize=12)
plt.yticks(LinRange(-0.001,0.003,5),fontsize=12);
plt.xticks([0,0.2,0.4,0.6],fontsize=12);




# ![RydbergBlockadeSubspace](../../../assets/RydbergBlockadeSubspace.png)

# Data for this evolution is shown above, where exact evolution in the full space is shown in black, and the truncated evolution in the subspace is shown by red dashed. It is clear that even though the Hilbert space is ``12\times`` smaller, the dynamics are faithfully reproduced, up to high frequency oscillations (inset) from adjacent atoms in the Rydberg state, similar to the high frequency oscillations of the 2 atom conditional blockade example. However, at longer times this subspace approximation fails to reproduce the full space (shown by divergence between black and red dashed), as the perturbative effects become relevant over longer timescales.

