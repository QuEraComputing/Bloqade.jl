# [Working with Subspace](@id subspace)


Due to the strong Rydberg interactions, only one Rydberg excitation is allowed within the blockade radius (see [Rydberg Blockade](@ref)). This is the called blockade constraint. 
In Bloqade, we take advantage of this effect by allowing users to run emulation in a truncated subspace.  This is done by throwing out states that violate the blockade constraint. 
This process could help us accelerate the emulation and reach a bigger system size. In this section, we will show how to create a blockade subsapce, create register in subspace, 
obtain Hamiltonian matrix in subspace, and run emulation in subspace. 



## Create Blockade Subspace

One can create a blockade subspace via `blockade_subspace` method if we know the atomic positions 

```@docs
blockade_subspace
```

For example, we can construct a blockade subspace of a square lattice
using the code below

```@example subspace
using Bloqade
atoms = generate_sites(SquareLattice(), 3, 3, scale=5.1)
space = blockade_subspace(atoms, 5.2)
```
where we have created a ``3*3`` square lattice with nearest neighbour atoms seperated by ``5.1 \mu m``. Then we have created a
blockade subpace with subspace radius being ``5.2 \mu m``. This means that if two atoms have a distance that is smaller than (or equal to)
``5.2 \mu m``, the blockade subspace does not contain states where both of them being in Rydberg states. The number of allowed states 
is 63 in the above case, which is much smaller than the full Hilbert space 512. 


Here `space` is of type `Subspace`

```@docs
Subspace
```

Other than by using atomic poisitons and specifying the subspace radius, we can also using a graph to create a subspace. In this case, the subsapce 
corresponds to the space composed by the [independent set](https://en.wikipedia.org/wiki/Independent_set_(graph_theory)) of this graph. Bloqade has an explicit function for this, by using graph as 
an input, and gives the subspace as the output. Here is an example code

```@example subspace
using graph

g = SimpleGraph(5)
edge_set = [(1,2), (1, 4), (2, 5), (3, 4)]
for (i,j) in edge_set
    add_edge!(g, i, j)
end 
space = independent_set_subspace(g)
```


## Create register in subspace

One can create the register in subspace by feeding the `space` object instead of an integer for the common register interfaces, e.g

```@repl subspace
zero_state(space)
product_state(bit"000_000_001", space)
```

Or if you have an existing state stored as a subtype of `AbstractVector`, we can also create the register using
the constructor

```@repl subspace
state = rand(ComplexF64, length(space))
reg = SubspaceArrayReg(state, space)
```

## Obtain matrix of Hamiltonian in subspace

The matrix projected in subspace of a given Hamiltonian can be obtained via
[`mat`](@ref) as well, e.g

```@repl subspace
h1 = rydberg_h(atoms; Δ=0.2, Ω=0.1)
mat(h1, space)
```


## Other Operations in Subspace

All other operations in subspace are the same as fullspace
case, e.g to run an emulation in subspace, one just need to use the
subspace register [`SubspaceArrayReg`](@ref) instead of the fullspace register [`ArrayReg`](@ref), the rest are all the same

```@example subspace
reg = zero_state(space)
prob = SchrodingerProblem(reg, 0.1, h1)
emulate!(prob)
statevec(reg)
```

The measurement in register with subspace is the same as that in the full space. 


## Create constrained local Hamiltonian using subspace 

Although we are able to emulate our Hamiltonian problem in the projected subspace, the long-range tail of  
the Rydberg interactions will be present in the subspace Hamiltonian. In certain cases, the user may want to 
cut off the long-range tail by only simulating a constrained short-range Hamiltonian, e.g. [the PXP model](https://arxiv.org/abs/2011.09486). 
In this case, we are able to using Bloqade to easily deal with such problem for arbitrary graph in arbitrary dimension. 

Let us take the PXP model in 1D as an example. We first create a 1D chain, and then generate a subspace by projecting out states without nearest-neighbour interactions. 

```@example subspace
using Bloqade
atoms = generate_sites(ChainLattice(), 10, scale=5.1)
space = blockade_subspace(atoms, 5.2)
register = product_state(bit"0101010101", space)
h = XTerm(length(atoms), Ω=1.0) - NTerm(length(atoms), Δ=0) 
prob = SchrodingerProblem(reg, 0.1, h)
emulate!(prob)
```
After creating the subspace, we have built a Hamiltonian by explicitly summing up the Rabi frequency term  and the detuning term by using [`XTerm`]@(ref) [`NTerm`]@(ref) respectively. 
In this way, we have created a local constraint Hamiltonian (without long-range tailed interactions). Futhermore, if we want to emulate 
quantum dynamics under this Hamiltonian, we just need to create a subspace register and emulate the system under the created Hamiltonian. 


