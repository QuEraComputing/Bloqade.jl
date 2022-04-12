# Registers

Bloqade follows the register interface in [Yao]. It uses register to 
represent a device and its internal quantum state.
As for simulators, the most used register types are `ArrayReg`
and `SubspaceArrayReg`. They both uses a dense array to store
the corresponding quantum state, except the later also stores
a subspace object.

## Basic Interfaces

The `ArrayReg` comes from `Yao`, and `SubspaceArrayReg` is a 
special register made for simulations in Rydberg blockade subspace.
You can use `zero_state` to create a register with its internal state to be zero state ``|00...00\ranlge``.

```@repl registers
using Bloqade
zero_state(5) # creates a 5-qubit register
```

to create more general product state, you can use `product_state`
function

```@repl registers
product_state(bit"10011")
```

where `bit"10011` is a special Julia string literal defined for
bitstrings.

## Operations

You can perform various operations on registers via standard Yao 
register interface, this includes applying operators using `apply!`, measure on the interal quantum state using `measure!` or
calculate an expectation
using `expect`.

