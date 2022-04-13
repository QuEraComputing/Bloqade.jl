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
