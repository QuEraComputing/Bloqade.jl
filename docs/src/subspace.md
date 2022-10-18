# [Working with Subspace](@id subspace)

Due to the strong Rydberg interactions, only one Rydberg excitation is allowed if the atoms are close to each other. 
We typically take this as the blockade radius, ``R_b``, which is the 
distance for which the Rydberg interaction is the same as the Rabi frequency, ``\Omega`` (see [Rydberg Blockade](@ref blockade)). This is the so-called blockade constraint. 

In Bloqade, we can take advantage of this effect by allowing users to run emulation in a truncated subspace, i.e., by throwing out states that violate the blockade constraint. 
This can help accelerate the simulation and enables simulation for a larger system size. 
In this section, we will show how to create a blockade subspace, create registers in the subspace, 
obtain the Hamiltonian matrix in the subspace, and run emulation in the subspace.

!!! note
    Note that the blockade radius ``R_b`` is the distance for which the Rydberg interaction is the same as the Rabi frequency, ``\Omega``. 
    For accurate simulation, however, it's not recommended to throw away the states that's close to the blockade radius. In other words, it's safer to set the subspace radius ``R_s`` to be smaller than ``R_b``, where we throw away the blockade violated states when the atoms are within ``R_s``. 
    For example, if we set ``R_s = 1/2 * R_b``, we will be throwing away states that have interaction energies at least ``2^6*\Omega``, which will be a good approximation. 
    See the [Rydberg Blockade](@ref blockade) page for recommendations on how to set ``R_b``, ``R_s``, and the atom lattice separation, ``a``.


## Create the Blockade Subspace

One can create a blockade subspace via the `blockade_subspace` method if we know the atomic positions: 

```@docs
blockade_subspace
```

For example, we can construct a blockade subspace of a square lattice
using the code below:

```@example subspace
using Bloqade
atoms = generate_sites(SquareLattice(), 3, 3, scale=5.1)
space = blockade_subspace(atoms, 5.2)
```
We first created a ``3*3`` square lattice with nearest neighbor atoms separated by ``5.1`` μm. 
Then we have created a
blockade subspace with the subspace radius, ``R_s``, being ``5.2`` μm. 
This means that if two atoms have a separation distance that is smaller than (or equal to)
``5.2`` μm, 
then the blockade subspace does not contain states where both of them being in the Rydberg states.
For the dictionary shown, the left is the new index of the states in the blockade subspace; 
in this case, there are 63 allowed states, which is much smaller than the full Hilbert space size 512.
The vectors on the right correspond to the base-10 representations of the states in bitstrings. 

Here `space` is of the type `Subspace`:

```@docs
Subspace
```

Other than using atomic positions and the subspace radius, we can also use a graph to create a subspace. In this case, the subspace 
corresponds to the space composed by the [independent sets](https://en.wikipedia.org/wiki/Independent_set_(graph_theory)) of this graph. Bloqade has an explicit function for this, by using a graph as 
an input, and produces the subspace as the output. Here is an example code:

```@example subspace
using Graphs

g = SimpleGraph(5)
edge_set = [(1,2), (1, 4), (2, 5), (3, 4)]
for (i,j) in edge_set
    add_edge!(g, i, j)
end 
space = independent_set_subspace(g)
```


## Create Registers in the Subspace

One can create a register in the subspace by feeding the `space` object instead of an integer for the common register interfaces, e.g.:

```@repl subspace
zero_state(space)
product_state(bit"000_000_001", space)
```

Alternatively, if you have an existing state stored as a subtype of `AbstractVector`, you can also create the register using
the constructor:

```@repl subspace
state = rand(ComplexF64, length(space))
reg = SubspaceArrayReg(state, space)
```

## Obtain the Hamiltonian Matrix in the Subspace

The matrix projected in the subspace of a given Hamiltonian can be obtained via
[`mat`](@ref) as well, e.g.:

```@repl subspace
atoms = generate_sites(SquareLattice(), 3, 3, scale=5.1);
space = blockade_subspace(atoms, 5.2);
h1 = rydberg_h(atoms; Δ=2.0*2π, Ω=1.0*2π)
mat(h1, space)
```


## Other Operations in the Subspace

All other operations in the subspace are the same as the fullspace
case. 
For example, to run an emulation in the subspace, one just need to use the
subspace register [`SubspaceArrayReg`](@ref) instead of the fullspace register [`ArrayReg`](https://docs.yaoquantum.org/dev/man/registers.html#Array-Registers).
The rest of the code are the same:

```@example subspace
reg = zero_state(space)
prob = SchrodingerProblem(reg, 1.0, h1)
emulate!(prob)
statevec(reg)
```

Measurements on the subspace register is the same as that in the full space. 


## Create Constrained Local Hamiltonians in the Subspace

Although we are able to emulate our Hamiltonian problem in the projected subspace, the long-range tail of  
the Rydberg interactions will be present in the subspace Hamiltonian. In certain cases, you may not want the long-range tail by only simulating a constrained short-range Hamiltonian, e.g. [the PXP model](https://arxiv.org/abs/2011.09486).
In this case, we can use Bloqade to easily deal with such problems for an arbitrary graph in an arbitrary dimension. 

Let us take the PXP model in 1D as an example. We first create a 1D chain and then generate a subspace by projecting out states that have nearest-neighbor interactions. 

```@example subspace
atoms = generate_sites(ChainLattice(), 10, scale=5.1)
space = blockade_subspace(atoms, 5.2)
register = product_state(bit"0101010101", space)
h = 2π * 4.0 * SumOfX(length(atoms)) - 2π * 1.0 * SumOfN(length(atoms))
prob = SchrodingerProblem(register, 0.2, h)
emulate!(prob)
```
After creating the subspace, we have built a Hamiltonian by explicitly summing up the Rabi frequency term  and the detuning term by using [`SumOfX`](@ref) and [`SumOfN`](@ref) respectively. 
In this way, we have created a local constraint Hamiltonian (without the long-range interaction tail). Futhermore, if we want to emulate 
quantum dynamics under this Hamiltonian, we just need to create a subspace register and emulate the system under the created Hamiltonian.


