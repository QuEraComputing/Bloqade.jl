# Hamiltonians

The simplest way of creating a Hamiltonian is via the `rydberg_h`
function

```@docs
rydberg_h
```

## The Hamiltonian Expressions

Inside this package, a more general definition of hamiltonians
are supported as Symbolic expressions, this
gives user the flexiblity to define various different kind of
hamltonian by simply writing down the expression.

There are currently 4 terms supported: [`RydInteract`](@ref),
[`XTerm`](@ref), [`NTerm`](@ref), [`ZTerm`](@ref). The terms
can be added up to compose a new hamiltonian, e.g

```@repl hamiltonian
using EaRyd
h = XTerm(5, 1.0) + ZTerm(5, 1.0)
```


## Convert Hamiltonian Expression to Matrices

The Hamiltonian expressions can be converted to matrices
via type conversion, e.g we can convert the above hamiltonian
to a `SparseMatrixCSC`

```@repl hamiltonian
using SparseArrays
SparseMatrixCSC(h)
```
