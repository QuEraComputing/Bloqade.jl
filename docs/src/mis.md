# [Maximum Independent Set](@id mis)

[Rydberg Blockade](@ref blockade) is one of the most important properties of neutral-atom quantum computing based on Rydberg states. 
It naturally encodes the [independent set](https://en.wikipedia.org/wiki/Independent_set_(graph_theory)) constraint. 
More specifically, Rydberg blockade implies that two atoms cannot be both excited to the Rydberg state ``|r\rangle`` if they are close to each other, 
whereas independent set constraint means two vertices cannot be both in the independent set when they are connected by an edge.
Thus, one can consider atoms in the Rydberg state as vertices in an independent set. See the proposal in [H. Pichler, et al.](https://arxiv.org/pdf/1808.10816.pdf) for more details. 

In particular, one can use the ground state of the [Rydberg Hamiltonian](@ref Hamiltonians) to encode the [maximum independent set problem](https://en.wikipedia.org/wiki/Independent_set_(graph_theory)), 
which is to find the largest independent set of a given graph. 
For a particular subclass of geometric graphs, the so-called unit disk graphs, 
the Rydberg Hamiltonian can encode the solution without any overhead in the number of qubits. 
In fact, see an experimental demonstration of quantum optimization of maximum independent set up to 289 qubits in [S. Ebadi, et al.](https://arxiv.org/abs/2202.09372).

In Bloqade, we provide several functions to support the simulation of solving independent set problems on neutral-atom quantum computers.
We list them here in this section.

## The Maximum Independent Set Problem
In graph theory, an independent set is a set of vertices in a graph such that no two of which are adjacent.
The problem of finding maximum independent sets (MIS) is NP-hard, i.e. it is unlikely to be solved in a time polynomial to the problem size.
However, for a graph with a small to intermediate size, various solution space properties, including finding the MIS size and enumerating all MISs, can be computed using the package [`GenericTensorNetworks`](https://github.com/QuEraComputing/GenericTensorNetworks.jl); please also refer to the related manual pages [the independent set problem](https://queracomputing.github.io/GenericTensorNetworks.jl/dev/tutorials/IndependentSet/) and [the maximal independent set problem](https://queracomputing.github.io/GenericTensorNetworks.jl/dev/tutorials/MaximalIS/).

A tutorial on how to solve the MIS problem using Bloqade is detailed in the [MIS tutorial](@ref mis-tutorial) page.

In the following, we list the APIs in the module `BloqadeMIS`, many of which support the simulation of variational quantum algorithms for solving the MIS problem.

## References

```@autodocs
Modules = [BloqadeMIS]
```
