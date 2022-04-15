# Observables

## Yao Register Interface

Bloqade register are Yao registers, thus they support
all the interfaces Yao supports, this including

- `measure`: measure the register.
- `expect`: calculate the expectation on given observable operator.

the observable operator can be constructed in the same way of hamiltonian expression, e.g

```@repl observable
using Bloqade
r = rand_state(10) # create a register with random state
expect(SumOfN(nsites=10), r) # calculates rydberg density
```

one can also `measure` on the register

```@repl observable
measure(r; nshots=10)
```

## Convenient Wrappers

Bloqade also provides a few convenient wrappers on commonly used observables

```@docs
rydberg_density
rydberg_corr
```

## Create general observables using operator expression

Bloqade make use of Yao's block system to represent
operator expression, e.g one can construct the Rydberg
correlation operator as

```@repl observable
corr(n, i, j) = chain(n, put(i=>Op.n), put(j=>Op.n))
```

You can make up any kind of quantum operator in this way
and use it with the [`expect`](@ref) or [`measure`](@ref)
function.

And because the hamiltonian is also an operator expression,
thus it can be used as an observable too

```@repl observable
r = rand_state(5)
pos = [(i, ) for i in 1:5]
h = rydberg_h(pos; Ω=0.1)
expect(h, r)
```

## Reference

Here are some common operators re-exported from `Yao`.

```@docs
expect
measure
X
Y
Z
```
