# Observables

As we discussed in [Register](@ref), Bloqade follows Yao Register Interface. 
Thus we support all the interfaces of register that Yao supports, which includes 

- `measure`: measure bitstring of the register.
- `expect`: calculate the expectation value of a given observable operator.

The observable operator can be constructed in the same way as that of the Hamiltonian expression, e.g

```@repl observable
using Bloqade
r = rand_state(10) # create a register with random state
expect(SumOfN(nsites=10), r) # calculates rydberg density
```

To measure the bitstrings on the register, one can just use `measure` on that register

```@repl observable
measure(r; nshots=10)
```

## Convenient Wrappers

Bloqade also provides a few convenient wrappers on commonly used observables for Rydberg systems, including the Rydberg density and two-point correlation functions. 

```@docs
rydberg_density
```
For example, if we want to measure the Rydberg density at each site or at a specific site, we can use the code below

```@repl observable
reg = rand_state(10)
n_each = rydberg_density(reg)
n_2 = rydberg_density(reg, 2)
```
To access two-point correlation functions, we can use the [`rydberg_corr`](@ref) below. 

```@docs
rydberg_corr
```

This function will output an matrix which stores correlation function for each two sites 

```@repl observable
nn_corr = rydberg_corr(Op.n, reg)
```

It worth mentioning that besides `Op.n`, other single-site operators includes the pauli operators `X`, `Y` and `Z`. 


## Create general observables using operator expression

Bloqade make use of Yao's block system to represent
operator expression, e.g one can construct the Rydberg
correlation operator as

```@repl observable
corr(n, i, j) = chain(n, put(i=>Op.n), put(j=>Op.n))
```

You can make up any kind of quantum operator in this way
and use it with the [`expect`](@ref) or [`measure`](@ref)
function, e.g. 

```@repl observable
reg = rand_state(10)
corr_XY = chain(10, put(2=>Op.n), put(4=>Op.X))
expect(corr_XY, reg) 
```
And because the Hamiltonian is also an operator expression,
thus it can be used as an observable too

```@repl observable
r = rand_state(5)
pos = [(i, ) for i in 1:5]
h = rydberg_h(pos; Ω=0.1)
expect(h, r)
```

Please refer to the [Hamiltonians](@ref) to see other operators that are supported to build the Haimiltonian. 

## Reference

Here are some common operators re-exported from `Yao`.

```@docs
expect
measure
X
Y
Z
```
