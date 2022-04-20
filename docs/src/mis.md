# [Maximum Independent Set](@id mis)

One of the most important property of Rydberg system is the Rydberg blockade,
which naturally embeds the maximum independent set problem into its ground state. In Bloqade, we provide several functions to work with such property and we will discuss them in this section.

## The maximum independent set problem
In graph theory, an independent set is a set of vertices in a graph, no two of which are adjacent.
The problem of finding maximum independent sets (MIS) is NP-hard, i.e. unlikely to be solved in a time polynomial to the problem size.
Even aproximating the MIS size ``\alpha(G)`` for a graph ``G=(V,E)`` is hard.
In this tutorial we study the MIS problem defined on diagonal-coupled unit-disk grid graphs (DUGG)(see [arxiv:2202.09372](https://arxiv.org/abs/2202.09372)).
Although these graphs have highly constraint topology, finding its MISs is NP-hard.
We show how to map the MIS problem on this graph to a Rydberg atom array hamiltonian,
and use two quantum algorithms, the standard QAOA and a variational quantum algorithm with specially parametrized waveform, to find maximum independent sets.

## References

```@autodocs
Modules = [BloqadeMIS]
```
