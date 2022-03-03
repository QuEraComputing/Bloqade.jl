# Rydberg Blockade

The Rydberg blockade mechanism is a crucial component of the operation of neutral atom computers. By including a term in the effective Hamiltonian which adds an interaction between adjacent Rydberg atoms, we may construct gates and nontrivial dynamics which create entanglement and correlation across the system.

The term is a van der Waals interaction where there is an energy shift if two adjacent atoms are in the Rydberg state. The Van der Waals interaction is

```
V_{ij} = \frac{C_6}{|r_i - r_j|^6}\hat n_i \hat n_j
```

where $\hat n_i=|r_i\rangle\langle r_i|$ is the number operator on the $i$th site, which is 1 if the atom is in the Rydberg state, and zero if the atom is in the ground state. The coefficeint $C_6 = 5393398 MHz/\mu m^6$ is the interaction strength; characteristically, this interaction has a strength $\approx 10MHz$ for two atoms seperated by $10\mu m$, a similar scale to the Rabi coupling between the ground and Rydberg state. Crucially, this can be seen as an energy shift on atom $j$, conditional on the state of atom $i$, and so can be used, in a loose sense, as a conditional logical sense. This is because the adjacent atom forces the atom to be off-resonant with the laser field.

This conditional drive can be seen given the following dynamics. Suppose two atoms are close to each other ($< 10 \mu m$) and so interact under Van der Waals. The left atom is either in a Rydberg state, or in the ground state, and the right atom is originally in the ground state. Then, a Rabi drive is applied to the right atom, which couples the atom's ground state to the Rydberg state. For this example, we choose a Rabi drive of $\Omega=4$MHz and distance between atoms $|r_i-r_j| = 7\mu m$, which gives a conditional detuning of $\approx 50$MHz. When the left atom is in the ground state (black, top), there are no interactions and the state of the right atom oscillates between the ground state and Rydberg state; for a particular choice of timing this executes a $\pi$ pulse, flipping the right atom from the ground to Rydberg state. However, when the left atom is in the Rydberg state (red, bottom), there is a large detuning on the right atom, which causes the transfer to the Rydberg state to be strongly suppressed. In this case, the right atom (up to perturbative corrections) is never in the Rydberg state.

![RydbergBlockade](https://user-images.githubusercontent.com/20091330/156390383-f6a2e81c-72f6-4cfe-baf2-ac052c89717e.png)

This conditional energy shift is the basis of the **Rydberg Blockade**. Because of the large energy shift from having two adjacent atoms in the Rydberg state, evolution from an atomic ground state with local Rabi couplings between ground and Rydberg is restricted to a low energy subspace of states where no two adjacent atoms are in the Rydberg state. Furthermore, because the interaction strength with distance is a large power law, one can define a characteristic scale set by the Rabi coupling. If two atoms are close such that the conditional detuning is much larger than the Rabi coupling, one can consider the atoms to be blockading each other, and both atoms cannot simultaniously be in the Rydberg state. In contrast, if two atoms are far away, the two atoms never blockade each other and both atoms can simultaniously be in the Rydberg state.

The allowed states are then **independent sets** of a **unit disk graph** defined by the positions of the atoms. A unit disk graph is a set of vertices and edges, where vertices represent every atom, and there are edges if the distance between vertices is less than some radius $|r_i-r_j|<R, given by the characteristic Rabi frequency and van der Walls detuning. The blockade constraint is incoded by independent sets, which are subsets  $S$ of vertices of the unit disk graph such that no two vertices in $S$ are adjacent. The number of independent sets in a graph is much smaller than the number of subsets of the graph, so correspondingly the Hilbert space of quantum evolution is much smaller. This subspace constraint can thus speed quantum simulation considerably, at very little cost.


![RydbergBlockadeSubspace](https://user-images.githubusercontent.com/20091330/156587739-6ed599d7-9223-4998-bcdf-4f47b5b9a717.png)

To emphisize the effectiveness of this independent set subspace, some example nonequilibrium dynamics are shown above, for a ring of 12 atoms seperated by $7\mu m$. The system is initialized into the ground state of all atoms, and driven by a $4MHz$ Rabi drive, which couples each atom's ground and Rydberg state. If the atoms were far apart and non-interacting, each atom would oscillate completely between its ground state and Rydberg state with a period of $4\pi/4 \mu$sec. However, because adjacent atoms shift to the Rydberg state concurrently, they are dynamically blockaded, causing the maximum Rydberg density to only be 1/2, corresponding to an antiferromagnetic $Z_2$ state.

Because of the blockade, we may choose the independent subspace as the Hilbert space of the system-- a reduction from $D = 2^N=4092$ to $D=322$, as shown by red dashed. It is clear that even though the Hilbert space is $12\times$ smaller, the dynamics are faithfully reproduced, up to high frequency oscillations (inset) from adjacent atoms in the Rydberg state, similar to the high frequency oscillations of the 2 atom conditional blockade example. However, at longer times this subspace approximation fails to reproduce the full space (shown by divergence between black and red dashed). Note that for this example, the distance between atoms was chosen to be in an intermediate regime (eg, at the edge of the unit disk), which reduces the blockade effect and amplifies the approximate nature of the blockade. If the atoms were chosen to be closer together (say, $5\mu m$) or the Rabi strength was reduced, the blockade approximation becomes much stronger.





