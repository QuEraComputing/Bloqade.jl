```@meta
CurrentModule = Bloqade
```

# [Emulation of the Quantum Dynamics](@id emulation)

After we create the Rydberg Hamiltonian and the register for storing the quantum information, we can 
simulate the quantum many-body dynamics. The coherent dynamics of
the system is governed by the Schrödinger equation. The emulation interface of Bloqade is designed in a  
define-and-run style.  With Bloqade, we have two major types of emulation:

- ODE-solver based emulation for most of the problems.
- Krylov based emulation for piecewise constant problems or QAOA-like problem.



### Define the ODE Emulation Problem

ODE solvers are the major backend we use for most of the exact quantum 
dynamics simulation. The ODE solvers for Bloqade are powered by the [DiffEq.jl package](https://diffeq.sciml.ai/).

Bloqade provides a special problem type [`SchrodingerProblem`](@ref)
that supports most of the 
[integrator interface](https://diffeq.sciml.ai/stable/basics/integrator/)
of `DiffEq`, and most of
[the solver options](https://diffeq.sciml.ai/stable/basics/common_solver_opts/). Here, we introduce some common use cases of the integrator and solver options. For more advanced usage of the solvers,
please refer to the above link.

```@docs
BloqadeODE.SchrodingerProblem
```

## Run ODE-based Emulation

To run the emulation, you need to define the exact evolution and solver
you would like to run with via [`BloqadeODE.SchrodingerProblem`](@ref), and then feed the corresponding object to the
`emulate!` function:

```@docs
emulate!
```

For example, we can simulate the quantum dynamics of a time-dependent Hamiltonian by the following codes:

```@example evolution
using Bloqade
atoms = generate_sites(SquareLattice(), 3, 3; scale=5.1);
clocks = [0.0, 0.1, 0.2, 0.3, 0.4];
wf = piecewise_constant(;clocks, values=2π*[1.0, 2.0, 3.0, 4.0]);
h = rydberg_h(atoms; Δ=2π*2.0, Ω=wf); # create the Hamiltonian 
reg = zero_state(length(atoms)); # create fullspace register
ev = SchrodingerProblem(reg, 0.3, h) # the second input is the total time
emulate!(ev)
```
With `emulate!`, the quantum state stored in `reg` has been updated to the state after the time evolution. 

If you want to do operations during the real-time evolution,
such as measuring observables, you can instead
using the integrator interface with `for` loop and with `TimeChoiceIterator` on your desired clocks, e.g.:

```@example evolution
integrator = init(ev, Vern8())
for _ in TimeChoiceIterator(integrator, [0.1, 0.25])
    ev.reg # state at selected time
    @show measure(ev.reg)[] # measure the state at each time
end
```

You can use any function on the `reg` object.  For calculating observables, 
please see the [Registers and Observables](@ref observables) section.

!!! tip

    Remember to make sure your operation does not mutate your state so that
    it won't affect the evolution itself, since the entire time evolution
    is simulated by mutating the state vector stored in
    `reg`. Thus, do not use any function that has a `!` in its name
    on the register `info.reg` unless you are certain about what you are
    doing.



## Choose an ODE Solver

One of the most powerful tool of the Julia ecosystem is the DiffEq ecosystem
that implements many different solvers. These solvers have different advantages and trade-offs. Since simulating a quantum many-body Schrödinger equation has some
special properties compared to a general ODE problem, we will discuss some
general heuristics in this section on how to choose a good ODE solver and
how to check if your simulation converges. Because many-body Schrödinger equation's stiffness is unknown, we will not be using stiff problem solvers, but instead using non-stiff problem algorithms or auto-switching algorithms.

For most of the cases, one can use the `VCABM` solver for a large system simulation. However, this method requires more memory, which can be a bottleneck with GPUs.

The `Vern` family is another set of solvers that is good for many-body
Schrödinger equation, such as `Vern6`, `Vern7`, and `Vern8`. They also
have relatively good memory usage when utilize GPUs.

For a more detailed list of solvers, please refer to [DiffEq:Full list of solvers](https://diffeq.sciml.ai/stable/solvers/ode_solve/#Full-List-of-Methods).
For more detailed explanation on ODE solvers, please refer to [DiffEq:Recommended Methods](https://diffeq.sciml.ai/stable/solvers/ode_solve/#Recommended-Methods).

If you are familiar with MATLAB or Python, you may wish to compare the same
methods that you use in MATLAB or Python; you can find the corresponding
solvers in Julia in [DiffEq:Translation from MATLAB/Python/R](https://diffeq.sciml.ai/stable/solvers/ode_solve/#Translations-from-MATLAB/Python/R).


## Adaptive Steps in ODE Solvers

Our ODE solvers use adaptive steps by default. It provides a significant speedup
compared to standard fixed-step methods (see [our benchmark here](https://github.com/yardstiq/bloqade_benchmarks/tree/main#results)).
However, if one expects to retrieve the results during the time evolution, e.g.,
plotting Rydberg densities with the evolution time, fixed-step methods are sometimes 
preferred.

More specifically, when the adaptive steps are turned on, the time steps might be large,
but if one is interested in measuring some observables in smaller time steps, then the adaptive step 
method will not produce accurate results for the finer time step, but instead output results at the specific adaptive steps. 
In this situation, it's better to use fixed-step methods at the clocks where the observables are measured.

One can use the code below to turn off the adaptive steps when setting up the [`SchrodingerProblem`](@ref):

```@example evolution
atoms = generate_sites(SquareLattice(), 3, 3; scale=5.1);
h = rydberg_h(atoms; Δ=2π*2.0, Ω= 2π*1.0); # create the Hamiltonian 
reg = zero_state(length(atoms)); 
prob = SchrodingerProblem(reg, 3.0, h, adaptive = false, dt=1e-3);
```

Here, we've specified the fixed time step as `dt = 1e-3`.
If one only expects the final state of the evolution,
or the intervals between each chosen clock is much larger than the maximum
step size, then adaptive steps are preferred.

### Define the Krylov Emulation Problem

The Krylov-based method expects time-independent Hamiltonians. One can define such a time evolution via [`KrylovEvolution`](@ref) object.

```@docs
KrylovEvolution
```

## Run Krylov-based Emulation

We can run the Krylov-based emulation in a similar way using [`emulate!`](@ref):

```@repl evolution
emulate!(KrylovEvolution(reg, clocks, h))
```

However, as its name suggests, the Krylov-based emulation is not a standard ODE problem that DiffEq  supports. Thus, it does not support the ODE problem interface, but it's more like a gate-based interface. For example, the object `KrylovEvolution` is iterable:

```@example evolution
for (step, reg, duration) in KrylovEvolution(reg, clocks, h)
    @show step
    @show reg
    @show duration
    println("==========")
end
```

## Krylov vs ODE Solvers

The [`KrylovEvolution`](@ref) uses the Krylov subspace methods to simulate the
time evolution of time-independent operators ``\exp(i\Delta t_i H)``, where ``\Delta t_i`` is the duration of time-independent Hamiltonian ``H`` at time ``t``. This method is more efficient when the evolution itself is a discrete evolution, e.g. in QAOA and with
[`piecewise_constant`](@ref) waveforms. As for other cases, ODE solvers
are usually more efficient than [`KrylovEvolution`](@ref).
