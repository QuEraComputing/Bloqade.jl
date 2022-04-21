```@meta
CurrentModule = Bloqade
```

# [Emulation of Shordinger Equation](@id emulation)

The dynamics of Rydberg system is described by a Shordinger Equation
of a hamiltonian. The emulation interface of Bloqade is design as 
define-and-run style.

## Defining the Emulation

There are two major types of emulation:

- ODE solver based emulation for most of the problems.
- Krylov based emulation for piecewise constant problems or QAOA-like problem.

### ODE Emulation

ODE solver is the major backend we uses for most of the exact quantum 
dynamics simulation. In Bloqade, this is powered by the [DiffEq.jl package](https://diffeq.sciml.ai/).

Bloqade provdes a special problem type [`SchrodingerProblem`](@ref)
that support most of the 
[integrator interface](https://diffeq.sciml.ai/stable/basics/integrator/)
of `DiffEq`, and most of
[the solver options](https://diffeq.sciml.ai/stable/basics/common_solver_opts/). Here we will introduce common use cases of the integrator and solver options. For more advanced usage of the solver,
please refer to the above link.

```@docs
BloqadeODE.SchrodingerProblem
```

### Krylov Emulation

The Krylov-based method expects time independent hamiltonians, one can define such evolution via [`KrylovEvolution`](@ref) object.

```@docs
KrylovEvolution
```

## Run ODE-based Emulation

To run the emulation, you need to define the exact evolution and solver
you would like to run with via either [`KrylovEvolution`](@ref) or
[`BloqadeODE.SchrodingerProblem`](@ref), then feed the corresponding object to
`emulate!` function

```@docs
emulate!
```

For example, we can simulate a constant hamiltonian

```@repl evolution
using Bloqade
atoms = generate_sites(SquareLattice(), 3, 3; scale=5.1)
clocks = [0.0, 0.1, 0.2, 0.3]
wf = piecewise_constant(;clocks, values=[1.0, 2.0, 3.0, 4.0])
h = rydberg_h(atoms; Δ=2.0, Ω=wf)
reg = zero_state(length(atoms)) # create fullspace register
ev = SchrodingerProblem(reg, 0.3, h)
emulate!(ev)
```

Or if you would like to do some operation during the evolution,
such as measure observables during the evolution, you can instead
using the integrator interface with `for` loop with `TimeChoiceIterator`
on your desired clocks, e.g

```@example evolution
integrator = init(ev, Vern8())
for _ in TimeChoiceIterator(integrator, [0.1, 0.25])
    ev.reg # state at selected time
    @show measure(ev.reg)[] # measure the state at each time
end
```

You can use any function on the `reg` object, for calculating observables
please see the [Observables](@ref) section.

!!! tip

    remember to make sure your operation does not mutate your state so that
    this won't effect the evolution itself, since the entire time evolution
    is simulated by keep mutating the state vector stored in
    `reg` which means do not use any function that has a `!` in its name
    on the register `info.reg` unless you are certain about what you are
    doing.


## Run Krylov-based Emulation

We can run the Krylov-based emulation in a similar way using [`emulate!`](@ref)

```@repl evolution
emulate!(KrylovEvolution(reg, clocks, h))
```

However, as its name points out, the Krylov-based emulation is not a standard ODE problem that DiffEq  supports, thus it does not support the ODE problem interface, but a more gate-like interface, e.g the object `KrylovEvolution` is iterable

```@example evolution
for (step, reg, duration) in KrylovEvolution(reg, clocks, h)
    @show step
    @show reg
    @show duration
    println("==========")
end
```

## Choosing ODE solver

One of the most powerful aspect of Julia ecosystem is the DiffEq ecosystem
that implements lots of different solvers. These solvers have different trade-offs, and since simulating many-body Schrodinger equation has some
special properties comparing to a general ODE problem, we will discuss some
general heurestics in this section about how to choose a good ODE solver and
how to check if your emulation converges. And because many-body Schrodinger equation's stiffness is unknown, we will not be using stiff problem solver, but using non-stiff problem algorithm or auto-switching algorithm.

Most of the case one can use `VCABM` solver for large system simulation, however this method is requires more memory which can be a bottleneck when
utilizing GPUs.

The `Vern` family is another set of solvers that is good for many-body
Schrodinger equation, such as `Vern6`, `Vern7` and `Vern8`, they also
have relatively good memory usage when utilize GPUs.

For more detailed list of solvers please refer to [DiffEq:Full list of solvers](https://diffeq.sciml.ai/stable/solvers/ode_solve/#Full-List-of-Methods).
For more detailed explaination on ODE solvers please refer to [DiffEq:Recommended Methods](https://diffeq.sciml.ai/stable/solvers/ode_solve/#Recommended-Methods).

If you come from MATLAB or Python, you may expecting to compare the same
method that you use in MATLAB or Python, you can find the corresponding
solvers in Julia in [DiffEq:Translation from MATLAB/Python/R](https://diffeq.sciml.ai/stable/solvers/ode_solve/#Translations-from-MATLAB/Python/R).


## Adaptive Steps in ODE solver

Our ODE solver uses adaptive steps by default. It provides a huge speedup
comparing to standard fixed step methods (see [our benchmark here](#)).
However, if one expects to retreive results during the time evolution, e.g
plotting Rydberg density changes with the time, fixed step method should be
preferred otherwise the ODE solver will give constant results between each
step

(add @Johnason's plot on ring emulation here)

On the other hand, if one only expects the final state of the evolution,
or the interval between each chosen clock is much larger than maximum
step size, adaptive step is preferred.

## Krylov vs ODE solver

The [`KrylovEvolution`](@ref) uses Krylov subspace methods to simulate the
time evolution as discrete time-independent time evolution operators ``\\exp(i\Delta t_i H(t))``, where ``\Delta t_i`` is the duration of time-independent Hamiltonian ``H(t)`` at time ``t``. This method is more efficient when the evolution itself is a discrete evolution, e.g QAOA,
[`piecewise_constant`](@ref) waveform. As for other cases, ODE solvers
are usually more efficient than [`KrylovEvolution`](@ref).
