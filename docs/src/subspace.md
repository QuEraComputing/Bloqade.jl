# [Working with Subspace](@id subspace)


Due to the strong Rydberg interactions, only one Rydberg excitation is allowed within the blockade radius (see [Rydberg Blockade](@ref)). This is the called blockade constraint. 
In Bloqade, we take advantage of this effect by allowing users to run emulation in a truncated subspace.  This is done by throwing out states that violate the blockade constraint. 
This process could help us accelerate the emulation and reach a bigger system size. In this section, we will show how to create a blockade subsapce, create register in subspace, 
obtain Hamiltonian matrix in subspace, and run emulation in subspace. 



## Create Blockade Subspace

One can create a blockade subspace via `blockade_subspace` method

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
blockade subpace with blockade radius being ``5.2 \mu m``. This means that if two atoms have a distance that is smaller than (or equal to)
``5.2 \mu m``, the blockade subspace does not contain states where both of them being in Rydberg states. The number of allowed states 
is 63 in the above case, which is much smaller than the full Hilbert space 512. 


Here `space` is of type `Subspace`

```@docs
Subspace
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

The matrix in subspace of a given Hamiltonian can be obtained via
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
