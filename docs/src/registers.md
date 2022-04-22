# Registers and Observables

Bloqade follows the register interface in [Yao](https://yaoquantum.org). It uses a register to 
represent a device and its internal quantum state.
For Bloqade, the most commonly used register types are `ArrayReg`
and `SubspaceArrayReg`. They both use a dense array to store
the corresponding quantum state. The only difference is that `SubspaceArrayReg` also stores
a subspace object. In this section, we will only cover how to create registers and perform operations on it in the full Hilbert space. For subspace operations, please refer to the [subspace](@ref) page. We will also discuss a few convenient wrappers on commonly used observables for Rydberg systems. 

In the analog mode, we use the states ``|g\rangle`` (ground state) and ``|r\rangle`` (Rydberg state) to encode a qubit. To be consistent with the standard language of qubits, we refer the states ``|g\rangle`` and ``|r\rangle`` as ``|0\rangle`` and ``|1\rangle`` here.

## Basic Interface

To create a register with all atoms to be in the ground state ``| 00..00 \rangle``, we can use 
the function [`zero_state`](@ref) by specifying the number of qubits

```@repl registers
using Bloqade
zero_state(5) # creates a 5-qubit register
```

To create a more general product state in the computational basis, one can use the [`product_state`](@ref) function by inputting its bitstring

```@repl registers
product_state(bit"10011")
```
where `bit"10011"` is a special Julia string literal defined for bitstrings.

One can also construct the [`ArrayReg`](@ref) or [`SubspaceArrayReg`](@ref) directly from arrays, e.g.

```@repl registers
ArrayReg(rand(ComplexF64, 2^5))
```

For a subspace register, one can create in the following way

```@repl registers
space = Subspace(5, [0, 2, 3, 7])
state = rand(ComplexF64, length(space))
reg = SubspaceArrayReg(state, space)
```

Here, ``[0, 2, 3, 7]`` are base-10 integer representations of the corresponding states with bitstrings.
For a more detailed guide on working with the subspace, please see
[subspace](@ref).

## Operations on registers

You can perform various operations on registers via the standard [Yao](https://yaoquantum.org)
register interface. This includes applying operators on quantum 
states by using [`apply!`](@ref), measuring bitstrings with a 
projection on the quantum state by using `measure!`, and
calculating the expectation value of certain observables by using 
[`expect`](@ref). To inspect the internal state of the register, one 
can use the [`statevec`](@ref) method


```@repl registers
reg = rand_state(3)
measure(reg; nshots=5)
expect(put(3,1=>X), reg)
statevec(reg)
apply!(reg, put(1=>X))
```

For more detailed introduction of the register interface, please
refer to [Yao:Array Registers](https://docs.yaoquantum.org/dev/man/array_registers.html) and [Yao:AbstractRegister](https://docs.yaoquantum.org/dev/man/registers.html).



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



## References

```@docs
arrayreg
apply!
measure!
expect
statevec
zero_state
rand_state
product_state
SubspaceArrayReg
set_zero_state!
X
Y
Z

```
