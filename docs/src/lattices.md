# Lattices

## Create a Lattice

Creating a lattice is very simple in EaRyd, e.g we can create a square lattice as following

```@repl quick-start
using EaRyd
generate_sites(SquareLattice(), 3, 3)
```

this generates the atom positions on a ``3\times 3`` square lattice using the [`generate_sites`](@ref)
function.

we support the following built-in lattice: [`SquareLattice`](@ref), [`KagomeLattice`](@ref), [`HoneycombLattice`](@ref), and more. Please refer to [Lattices](@ref) for more detailed guide of lattice related operation.

!!! tip
    The lattice in EaRyd is actually defined as general Bravis lattice. You can create
    your own by defining the corresponding lattice vector. Check [Bravais Lattice](@ref bravais-lattice)

## Bravais Lattices

## Plotting Lattices
