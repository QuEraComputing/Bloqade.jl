# [Maximum Independent Set](@id mis)

One of the most important property of Rydberg system is the Rydberg blockade,
which naturally embeds the [maximum independent set problem](https://en.wikipedia.org/wiki/Independent_set_(graph_theory)) into its ground state. In Bloqade, we provide several functions to work with such property and we will discuss them in this section.

## The maximum independent set problem
In graph theory, an independent set is a set of vertices in a graph, no two of which are adjacent.
The problem of finding maximum independent sets (MIS) is NP-hard, i.e. unlikely to be solved in a time polynomial to the problem size.
However, for a graph with a small to intermediate size, various solution space properties can be solved with the package [`GenericTensorNetworks`](https://github.com/QuEraComputing/GenericTensorNetworks.jl), related manual pages are [the independent set problem](https://queracomputing.github.io/GenericTensorNetworks.jl/dev/tutorials/IndependentSet/) and [the maximal independent set problem](https://queracomputing.github.io/GenericTensorNetworks.jl/dev/tutorials/MaximalIS/).

In the following, we list the APIs in module `BloqadeMIS`, many of them are for emulating the variational quantum algorithm for solving the MIS problem in manual page [The Maximum Independent Set Problem](@ref).

## References

```@autodocs
Modules = [BloqadeMIS]
```
