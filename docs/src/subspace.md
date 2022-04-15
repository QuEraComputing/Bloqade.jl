# [Working with Subspace](@id subspace)

## Blockade Subspace

The blockade subspace of Rydberg system can be used to reduce Hilbert space
size so that one can approximate large system using relatively small space.

One can create a blockade subspace via `blockade_subspace` method

```@docs
blockade_subspace
```

For example, we can construct a blockade subspace of a square lattice
with random dropout.

```@example subspace
using Bloqade
atoms = generate_sites(SquareLattice(), 3, 3, scale=5.1)
space = blockade_subspace(atoms, 5.1)
```

here `space` is of type `Subspace`

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

## Run emulation in subspace

To run an emulation in subspace, one just need to use the
subspace register [`SubspaceArrayReg`](@ref) instead of the fullspace register [`ArrayReg`](@ref), e.g

```@example subspace
reg = zero_state(space)
prob = SchrodingerProblem(reg, 0.1, h1)
emulate!(prob)
statevec(reg)
```
