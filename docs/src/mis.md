# [Maximum Independent Set](@id mis)

One of the most important property of Rydberg system is the Rydberg blockade,
which naturally embeds the maximum independent set problem into its ground state. In Bloqade, we provide several functions to work with such property and we will discuss them in this section.

## Blockade Subspace

The blockade subspace of Rydberg system can be used to reduce Hilbert space
size so that one can approximate large system using relatively small space.

One can create a blockade subspace via `blockade_subspace` method

```@docs
blockade_subspace
```

## References

```@autodocs
Modules = [BloqadeMIS]
```
