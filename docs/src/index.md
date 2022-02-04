```@meta
CurrentModule = EaRyd
```

# EaRyd

[![Build Status](https://github.com/Happy-Diode/EaRyd.jl/workflows/CI/badge.svg)](https://github.com/Happy-Diode/EaRyd.jl/actions)
[![Coverage Status](https://coveralls.io/repos/github/Happy-Diode/EaRyd.jl/badge.svg?branch=master&t=p1FNvJ)](https://coveralls.io/github/Happy-Diode/EaRyd.jl?branch=master)

Welcome to the QuEra **E**mul**a**tor for **Ryd**berg System documentation page!

## Installation

```@raw html
<p>
EaRyd is a &nbsp;
    <a href="https://julialang.org">
        <img src="https://raw.githubusercontent.com/JuliaLang/julia-logo-graphics/master/images/julia.ico" width="16em">
        Julia Language
    </a>
    &nbsp; package. To install EaRyd,
    please <a href="https://docs.julialang.org/en/v1/manual/getting-started/">open
    Julia's interactive session (known as REPL)</a> and press <kbd>]</kbd> key in the REPL to use the package mode, then type the following command
</p>
```

First add the QuEra Julia registry via

```julia
pkg> registry add https://github.com/Happy-Diode/Miskatonic.git
```

For stable release

```julia
pkg> add EaRyd
```

For current master

```julia
pkg> add EaRyd#master
```

## Rydberg System

Our Rydberg Emulator simulate the following interacting Rydberg Hamiltonian, 

```math
\frac{\mathcal{H}}{\hbar} = \sum_i \frac{\Omega_i}{2} \sigma_x^i - \sum_i \Delta_i n_i + \sum_{i < j} V_{ij} n_i n_j.
```

Here ``\Delta_i`` are the detunings of the driving lasers from the Rydberg state;  ``\sigma_x^i = |g_i\rangle \langle r_i| + |r_i\rangle \langle g_i|`` describes the coupling between the ground state ``|g_i\rangle`` and the Rydberg state ``|r_i\rangle`` of an atom at position ``i``, driven at Rabi frequencey ``\Omega_i``;  ``n_i = |r_i\rangle \langle r_i|``, and ``V_{ij} = C/r_{ij}^6`` is the interactions between atom ``i`` and ``j``, where ``C`` is the interacting constant that depends on particular Rydberg atoms, and ``\hbar`` is the reduced Plank's constant. 

## Run a Simple Emulation of Rydberg System

Here we will show a simple example about simulating quantum many-body dynamics governed by such a Hamiltonian. 

We start by loading the Emulator Module

```@repl quick-start
using EaRyd
```

As one can see from the Rydberg Hamiltonian, the interactions between Rydberg atoms depend on their positions. EaRyd provides several built-in lattice structures for specifying the atom positions. For instance, we can use the following codes to quickly generate a chain of 10 atoms in 1D 

```@repl quick-start
nsites = 10;
atoms = generate_sites(ChainLattice(), nsites, scale=5.72)
```
We have set the distance between nearest neighbor atoms to be 5.72 ``\mu m``. Note that the default unit of length in EaRyd is ``\mu m``. Correspondingly, the default value for interacting constant ``C`` is ``2π * 858386`` MHz ``\times \mu m^6`` for Rydberg atom ``^{87}Rb`` and ``70 s`` Rydberg state. 

We will set both ``\Omega`` and ``\Delta`` to be a constant. The default energy unit (also for ``\Omega`` and ``\Delta``) is MHz in EaRyd. Since all the parameters specified, we can create a interacting Rydberg Hamiltonian by using [`rydberg_h`](@ref), 

```@repl quick-start
h = rydberg_h(atoms;C = 2π * 858386, Ω=4π, Δ=0)
```

We use the following code to create an initial state with all the atoms in the ground state by using [`zero_state`](@ref)

```@repl quick-start
init = zero_state(10)
```

We are interested in calculating the quench dynamics of Rydberg system starting with the above initial state and under the evolution of the defined Hamiltonian. Suppose we only want to get observable expectation values at the final time step, we can first create the problem and then directly emulate the problem.

```@repl quick-start
prob = ODEEvolution(init, 1.6, h)
emulate!(prob)
```
Here we have choosen the ODE solver [`ODEEvolution`](@ref) and set the total evolution time to be 1.6 ``\mu s`` (the default unit for time is ``\mu s`` in EaRyd). 

After emulating the problem, we can measure the Rydberg density at each site for the final state 

```@repl quick-start
densities = map(1:nsites) do i
    real(expect(put(nsites, i=>Op.n), prob.reg))
end
```
`prob.reg` is the register storing the final state after the time-evolution. 


## Looking for Help?

- check the slack channel [#julia](https://quera-workspace.slack.com/archives/C011C12GXRD)
- if not urgent, ask questions in [discussions](https://github.com/Happy-Diode/EaRyd.jl/discussions)

## Have Suggestions or Interested in Contribution?

- check the slack channel [#q-emulator](https://quera-workspace.slack.com/archives/C01MKUATZRD) for meetings and discussions
- [file an issue](https://github.com/Happy-Diode/EaRyd.jl/issues/new) to report a bug or request a feature
