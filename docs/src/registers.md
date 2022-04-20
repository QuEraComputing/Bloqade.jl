# Registers

Bloqade follows the register interface in [Yao](https://yaoquantum.org). It uses register to 
represent a device and its internal quantum state.
As for our Rydberg emulator, the most commonly used register types are `ArrayReg`
and `SubspaceArrayReg`. They both use a dense array to store
the corresponding quantum state. The only difference is that `SubspaceArrayReg` also stores
a subspace object. In this section, we will only cover how to create and operate registers in the full space. For subspace operations, please refer to [subspace](@ref). 

## Basic Interfaces


To create a register with its internal state to be Rydberg ground state ``| 00..00 \rangle``, we can simply use 
the function [`zero_state`](@ref) by specifying the number of qubits

```@repl registers
using Bloqade
zero_state(5) # creates a 5-qubit register
```

To create a more general polarized product state, you can use the [`product_state`](@ref) function by inputing its bit-string

```@repl registers
product_state(bit"10011")
```
where `bit"10011` is a special Julia string literal defined for bitstrings.

One can also construct the [`ArrayReg`](@ref) or [`SubspaceArrayReg`](@ref) directly from arrays, e.g

```@repl registers
ArrayReg(rand(ComplexF64, 1<<5))
```

Or for subspace register, one can create in the following way

```@repl registers
space = Subspace(5, [0, 2, 3, 7])
state = rand(ComplexF64, length(space))
reg = SubspaceArrayReg(state, space)
```

For more detailed guide about working with subspace, please see
[subspace](@ref).

## Operations

You can perform various operations on registers via standard [Yao](https://yaoquantum.org)
register interface. This includes applying operators on quantum 
states by using [`apply!`](@ref), measuring bitstrings or certain observables with 
projection on the quantum state by using `measure!`, and
calculating the expectation value of certain observables by using 
[`expect`](@ref). To inspect the internal state of the register, one 
can use the [`statevec`](@ref) method


```@repl registers
reg = rand_state(3)
measure(reg)
expect(put(1=>X), reg)
statevec(reg)
apply!(reg, put(1=>X))
```

For more detailed introduction of register interface, please
refer to [Yao:Array Registers](https://docs.yaoquantum.org/dev/man/array_registers.html) and [Yao:AbstractRegister](https://docs.yaoquantum.org/dev/man/registers.html).


## References

```@docs
arrayreg
apply!
measure!
statevec
zero_state
rand_state
product_state
SubspaceArrayReg
set_zero_state!
```
