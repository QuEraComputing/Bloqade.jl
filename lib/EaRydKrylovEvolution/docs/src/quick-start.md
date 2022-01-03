# Quick Start

```@setup example
using EaRydKrylovEvolution
```

Here are some quick example of doing some common tasks. First you need to
import some functions from this package via

```julia
using EaRydKrylovEvolution
```

## Create List of Atom Positions

When we need to emulate the interaction term, we will need to create a list of positions,
we support n-dimensional atom position via [`RydAtom`](@ref)

```@docs
RydAtom
```

A convenient function [`square_lattice`](@ref) is provided to create a list of atoms
on a square lattice

```@repl example
square_lattice(10, 0.8)
```

which creates a list of `10` atoms on square lattice with filling factor `0.8`.

!!! warn
    The `square_lattice` function should be considered temprorary convenient function
    it may be deprecated by more generic implementations in the future.

## Creating Rydberg Hamiltonians

Emulating the dynamics of a Hamiltonian is in principle very simple,
which is basically just about evaluating ``exp(-itH) \cdot state``,
the emulator has a symbolic representation of the hamiltonian to support
composing and simulating different kind of hamiltonian efficiently.

First, we can construct a hamiltonian only contains Pauli X term [`XTerm`](@ref).

```@docs
XTerm
```

We can create a 5-atom [`XTerm`](@ref) using the constructor, it has a global
constant rabi frequency `1.0`,

```@repl example
XTerm(5, 1.0)
```

we can also put a phase in it.

```@repl example
XTerm(5, 1.0, 1.2)
```

We can also create other hamiltonian terms such as [`RydInteract`](@ref), [`ZTerm`](@ref),
[`NTerm`](@ref) using their constructors, and then we can compose them together using `+`
or `-` operator

```@repl example
XTerm(5, 1.0, 1.2) + NTerm(5, 2.2)
XTerm(5, 1.0, 1.2) - NTerm(5, 2.2)
```

The parameters of the hamiltonian can also be a list of numbers that specifies the parameter
of each atom, in this case, you don't need to write the number of atoms explicitly

```@repl example
XTerm(rand(5))
```

The parameters can also be a function that specifies a pulse shape

```@repl example
h = XTerm(5, sin)
```

If the parameters are functions, we can use the hamiltonian term as a function of time `t`

```@repl example
h(1.2)
```

## Emulate Rydberg dynamics

Given a hamiltonian

```@repl example
atoms = square_lattice(5, 0.8);
h = RydInteract(atoms) + XTerm(length(atoms), 1.2) - NTerm(length(atoms), 1.1)
```

we can simulate its dynamics via [`emulate`](@ref)

```@repl example
r = emulate([2e-3], [h])
```

This will return us a register of type [`ArrayReg`](https://docs.yaoquantum.org/dev/man/array_registers.html), which is the type used in the general purpose quantum computing framework in Julia [Yao](https://yaoquantum.org).

Thus we can use all other functionalities implemented for `ArrayReg` from `Yao` directly, such as `measure`

```@repl example
using Yao
measure(r; nshots=10)
```

or we can calculate expectations of other operators

```@repl example
expect(put(5, 1=>X), r)
```

or other hamiltonians

```@repl example
using YaoExtensions
expect(heisenberg(5), r)
```

The more advanced usage of [`emulate`](@ref) function can be found as following

```@docs
emulate
```

The emulator also provides some convenient functions such as [`mean_rydberg`](@ref)

```@docs
mean_rydberg
```

this would return the mean size of vertices.

```@repl example
mean_rydberg(r)
```

## Emulate Rydberg dynamics in blockade subspace

In most cases, we can simulate the Rydberg dyanmics in blockade subspace,

```@repl example
atoms = square_lattice(10, 0.8)
h = RydInteract(atoms) + XTerm(length(atoms), 1.2) - NTerm(length(atoms), 1.1)
space = blockade_subspace(atoms, 1.0)
```

This would have a much smaller Hilbert space size comparing to ``2^{10}``

```@repl example
length(space)
```

Similarly we can create run the simulation in the subspace

```@repl example
r = emulate(space, [1e-3], [h])
```

This would return a `RydbergReg` that stores the amplitudes in the subspace,
but it supports the same interface as `ArrayReg`

```@repl example
measure(r; nshots=5)
```

## Solving Maximum Independent Set Problem via Adiabatic Schedule

TBD.

## Solving Maximum Independent Set Problem via QAOA

TBD.