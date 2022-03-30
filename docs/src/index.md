```@meta
CurrentModule = EaRyd
```

# EaRyd

[![Coverage Status](https://coveralls.io/repos/github/Happy-Diode/EaRyd.jl/badge.svg?branch=master&t=p1FNvJ)](https://coveralls.io/github/Happy-Diode/EaRyd.jl?branch=master)

Welcome to the documentation page for the QuEra **E**mul**a**tor for **Ryd**berg System.

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
    Julia's interactive session (known as REPL)</a>, press <kbd>]</kbd> key in the REPL to use the package mode, and then type the following command
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

## What does the Rydberg Emulator do?

This Rydberg Emulator package simulates the time evolution of a quantum state under the Schrödinger equation where the Hamiltonian is the interacting Rydberg Hamiltonian `` \mathcal{H} ``, 

```math
i \hbar \dfrac{\partial}{\partial t} | \psi \rangle = \mathcal{H} | \psi \rangle,  \\

\frac{\mathcal{H}}{\hbar} = \sum_j \frac{\Omega_j(t)}{2} \left( e^{i \phi_j(t) } | 1_j \rangle  \langle r_j | + e^{-i \phi_j(t) } | r_j \rangle  \langle 1_j | \right) - \sum_j \Delta_j(t) n_j + \sum_{j < k} V_{jk} n_j n_k.
```

Here ``\Omega_j``, ``\phi_j``, and ``\Delta_j``  denote the Rabi frequency, laser phase, and the detuning of the driving laser field on atom (qubit) ``j`` coupling the two states  ``| 1_j \rangle `` (ground state) and `` | r_j \rangle `` (Rydberg state). The number operator ``n_j = |r_j\rangle \langle r_j|``, and ``V_{jk} = C/|\overrightarrow{\mathbf{r_j}} - \overrightarrow{\mathbf{r_k}}|^6`` describes the Rydberg interaction between atoms ``j`` and ``k`` where ``\overrightarrow{\mathbf{r_j}}`` denotes the position of the atom ``j``; ``C`` is the Rydberg interaction constant that depends on the particular Rydberg state used. For the emulator, the default ``C = 2\pi \times 862690 \text{ MHz μm}^6`` for ``|r \rangle = 70S_{1/2}`` of the ``^{87}``Rb atoms. ``\hbar`` is the reduced Planck's constant.

Starting from the initial quantum state ``| \psi_{\text{ini}} \rangle``, the emulator can simulate the time evolution of the quantum state under the time-dependent Hamiltonian ``\mathcal{H}(t)``, given the atom positions, the time-dependent profiles for  ``\Omega_j``, ``\phi_j`` and ``\Delta_j``. The emulator can then produce the real-time evolution of quantum state ``| \psi(t) \rangle`` and one can then measure different observables on such a state.

The default units for various quantities are 

| Quantity      | Default Unit |
| :---:         |    :----:   |
| Length        |  μm         |
| Time          |  μs         |
| ``\Omega``    |  MHz        |
| ``\phi``      |  rad        |
| ``\Delta``    |  MHz        |

## Run a Simple Emulation of Rydberg System

Here, we show a simple example of simulating quantum many-body dynamics governed by the Rydberg Hamiltonian. 

We start by loading the Emulator Module

```@repl quick-start
using EaRyd
```

As one can see from the Rydberg Hamiltonian, the interactions between Rydberg atoms depend on their positions. EaRyd provides several built-in [Lattices](@ref) structures for specifying the atom positions. For instance, we can use the following codes to quickly generate a chain of 10 atoms in 1D: 

```@repl quick-start
nsites = 10;
atoms = generate_sites(ChainLattice(), nsites, scale = 5.74)
```
We have set the distance between nearest neighbor atoms to be 5.74 μm. Note that the default unit of length is μm as shown in the table above.

Let's set both ``\Omega`` and ``\Delta`` to be constants. Since all the parameters are specified, we can now create an interacting Rydberg Hamiltonian by using [`rydberg_h`](@ref), 

```@repl quick-start
h = rydberg_h(atoms; Ω = 2π * 4, Δ = 0)
```

For creating more complicated waveforms for ``\Omega`` and ``\Delta`` and the supported utilities, please refer to the [Waveforms](@ref) page.

Let's create an initial state with all the atoms in the ground state by using [`zero_state`](@ref)

```@repl quick-start
reg = zero_state(10)
```

We are interested in measuring observables of the final quantum state of the Rydberg system starting from the initial state and evolving under the Rydberg Hamiltonian over some time duration. We can first create the problem and then directly emulate the time evolution.

```@repl quick-start
prob = SchrodingerProblem(reg, 1.6, h)
integrator = init(prob, Vern8())
emulate!(prob);
```
Here we have chosen the ODE solver [`ODEEvolution`](@ref) and set the total evolution time to be 1.6 μs.

After simulating the time evolution and get the final state, we can measure the Rydberg population at each site for the final state 

```@repl quick-start
rydberg_populations = map(1:nsites) do i
    real(expect(put(nsites, i=>Op.n), prob.reg))
end
```
`prob.reg` is the register storing the final state after the time evolution.


## Looking for Help?

- check the slack channel [#julia](https://quera-workspace.slack.com/archives/C011C12GXRD)
- if not urgent, ask questions in [discussions](https://github.com/Happy-Diode/EaRyd.jl/discussions)

## Have Suggestions or Interested in Contribution?

- check the slack channel [#q-emulator](https://quera-workspace.slack.com/archives/C01MKUATZRD) for meetings and discussions
- [file an issue](https://github.com/Happy-Diode/EaRyd.jl/issues/new) to report a bug or request a feature
