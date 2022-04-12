```@meta
CurrentModule = Bloqade
```

# Emulation of Shordinger Equation

The dynamics of Rydberg system is described by a Shordinger Equation
of a hamiltonian. There are two methods for solving this equation,
one is via Krylov subspace projection for time independent hamiltonian,
the other is ODE-basd solver.

## Krylov Based Solver

The Krylov-based method expects time independent hamiltonians, one can define such evolution via [`KrylovEvolution`](@ref) object.

```@docs
KrylovEvolution
```

!!! tip

    one can use [`trotterize`](@ref) function to discretize a time-dependent
    hamiltonian into a list of time indepdent hamiltonian to use 
    Krylov-based solver for time-dependent problem if preferred.

## ODE Based Solver

The ODE-based method expects time-dependent hamiltonians, one can define such evolution via [`ODEEvolution`](@ref) object.

```@docs
ODEEvolution
```

## Run Emulation

To run the emulation, you need to define the exact evolution and solver
you would like to run with via either [`KrylovEvolution`](@ref) or
[`ODEEvolution`](@ref), then feed the corresponding object to
`emulate!` function

```@docs
emulate!
```

For example, we can simulate a constant hamiltonian

```@repl evolution
using Bloqade
atoms = generate_sites(SquareLattice(), 3, 3; scale=5.1)
ds = rand(3) # durations
hs = [rydberg_h(atoms;Δ=2.0, Ω=1.0) for _ in 1:3]
reg = zero_state(length(atoms)) # create fullspace register
ev = KrylovEvolution(reg, ds, hs)
emulate!(ev)
```

Or if you would like to do some operation during the evolution,
such as measure observables during the evolution, you can instead
write the `for` loop

```@example evolution
ev = KrylovEvolution(reg, ds, hs)
for info in ev
    @info "running emulation" step=info.step duration=info.duration h=info.hamiltonian reg=info.reg
end
```

You can use any function on the `reg` object.

!!! tip

    remember to make sure your operation does not mutate your state so that
    this won't effect the evolution itself, since the entire time evolution
    is simulated by keep mutating the state vector stored in
    `reg` which means do not use any function that has a `!` in its name
    on the register `info.reg` unless you are certain about what you are
    doing.
