```@meta
CurrentModule = Bloqade
```

# Bloqade

```@raw html
<p>
Welcome to the documentation page for Bloqade, a &nbsp;
    <a href="https://julialang.org">
        <img src="https://raw.githubusercontent.com/JuliaLang/julia-logo-graphics/master/images/julia.ico" width="16em">
        Julia Language
    </a>
    &nbsp; package for the simulation of quantum computation and quantum dynamics based on neutral-atom architectures.
```

Neutral-atom quantum computers have two major modes of computation: the first mode is a "digital mode" to do universal, digital quantum computation that uses two ground states ``|0\rangle`` and ``|1\rangle`` to encode the qubit, which has long coherence time, and one Rydberg state ``|r\rangle`` to entangle the qubits; the second mode is an "analog mode" as a programmable quantum simulator that uses one ground state ``|1\rangle`` and one Rydberg state ``|r\rangle``, where the quantum dynamics is governed by a Rydberg Hamiltonian ``\hat{\mathcal{H}}``.

Currently, Bloqade enables the easy design and fast execution of quantum dynamics in analog mode,  based on the neutral-atom quantum computing architecture. Besides fast full Hilbert-space simulation on CPUs, the main features include the design of arbitrary-layout quantum registers ([Lattices](@ref)), easy waveform generation ([Waveforms](@ref)), the simulation in subspace constrained by the Rydberg blockade ([subspace](@ref)), faster GPU-accelerated simulation ([CUDA Acceleration](@ref)), and more.

## Installation

```@raw html
<p>
To install Bloqade,
    please <a href="https://docs.julialang.org/en/v1/manual/getting-started/">open
    Julia's interactive session (known as REPL)</a>, press <kbd>]</kbd> key in the REPL to use the package mode, and then
</p>
```

add the QuEra Julia registry via

```julia
pkg> registry add https://github.com/Happy-Diode/Miskatonic.git
```

For the stable release, type

```julia
pkg> add Bloqade
```

Or for the current master,

```julia
pkg> add Bloqade#master
```

For a more advanced installation guide, please see the [Installation](@ref install) page.

## What does Bloqade do?


In the analog mode, Bloqade simulates the time evolution of a quantum state under the Schrödinger equation where the Hamiltonian is the interacting Rydberg Hamiltonian ``\hat{\mathcal{H}}``, 

```math
i \hbar \dfrac{\partial}{\partial t} | \psi \rangle = \hat{\mathcal{H}}(t) | \psi \rangle,  \\

\frac{\mathcal{H}(t)}{\hbar} = \sum_j \frac{\Omega_j(t)}{2} \left( e^{i \phi_j(t) } | 1_j \rangle  \langle r_j | + e^{-i \phi_j(t) } | r_j \rangle  \langle 1_j | \right) - \sum_j \Delta_j(t) \hat{n}_j + \sum_{j < k} V_{jk} \hat{n}_j \hat{n}_k.
```

Following the atomic physics nomenclature, ``\Omega_j``, ``\phi_j``, and ``\Delta_j``  denote the Rabi frequency, laser phase, and the detuning of the driving laser field on atom (qubit) ``j`` coupling the two states  ``| 1_j \rangle `` (ground state) and `` | r_j \rangle `` (Rydberg state); ``\hat{n}_j = |r_j\rangle \langle r_j|`` is the number operator, and ``V_{jk} = C_6/|\overrightarrow{\mathbf{r_j}} - \overrightarrow{\mathbf{r_k}}|^6`` describes the Rydberg interaction (van der Waals interaction) between atoms ``j`` and ``k`` where ``\overrightarrow{\mathbf{r_j}}`` denotes the position of the atom ``j``; ``C_6`` is the Rydberg interaction constant that depends on the particular Rydberg state used. For Bloqade, the default ``C_6 = 862690 \times 2\pi \text{ MHz μm}^6`` for ``|r \rangle = \lvert 70S_{1/2} \rangle`` of the ``^{87}``Rb atoms; ``\hbar`` is the reduced Planck's constant.

Starting from an initial quantum state ``| \psi_{\text{ini}} \rangle``, Bloqade simulates its time evolution under the Hamiltonian ``\hat{\mathcal{H}}(t)``, given the qubit positions and the time-dependent profiles for  ``\Omega_j``, ``\phi_j``, and ``\Delta_j``. Bloqade then outputs the real-time-evolved state ``| \psi(t) \rangle``, which can then be used for measuring different observables.

The default units for various quantities are 

| Quantity      | Default Unit |
| :---:         |    :----:   |
| Length        |  μm         |
| Time          |  μs         |
| ``\Omega``    |  MHz        |
| ``\phi``      |  rad        |
| ``\Delta``    |  MHz        |


## First steps

Let's try a simple example of simulating quantum many-body dynamics governed by the Rydberg Hamiltonian. 

We start by loading the Bloqade Module

```@repl quick-start
using Bloqade
```

As one can see from the Rydberg Hamiltonian, the interactions between Rydberg atoms depend on their positions. Bloqade provides several built-in [Lattices](@ref) structures for specifying the atom positions. For instance, we can use the following codes to quickly generate a chain of 10 atoms in 1D: 

```@repl quick-start
nsites = 10;
atoms = generate_sites(ChainLattice(), nsites, scale = 5.74)
```
We have set the distance between nearest-neighbor atoms to be 5.74 μm. Note that the default unit of length is μm as shown in the table above.

Let's set both ``\Omega`` and ``\Delta`` to be constants (and ``\phi = 0``). Since all the variable parameters in the Hamiltonian are specified, we can now create an interacting Rydberg Hamiltonian by using [`rydberg_h`](@ref), 

```@repl quick-start
h = rydberg_h(atoms; Ω = 4 * 2π, Δ = 0)
```

To create more complicated waveforms for ``\Omega`` and ``\Delta`` and find the supported utilities, please refer to the [Waveforms](@ref) page.

Let's create an initial state with all the atoms in the ground state by using [`zero_state`](@ref)

```@repl quick-start
reg = zero_state(10)
```

We are interested in measuring observables of the final quantum state of the Rydberg system starting from the initial state and evolving under the Rydberg Hamiltonian over some time duration. We can first create the problem and then directly simulate the time evolution.

```@repl quick-start
prob = SchrodingerProblem(reg, 1.6, h)
integrator = init(prob, Vern8());
emulate!(prob);
```
Here, we have chosen the ODE-based solver (`Vern8()`) by using [`SchrodingerProblem`](@ref) and set the total evolution time to be 1.6 μs.

After simulating the time evolution and get the final state, we can measure the Rydberg population at each site for the final state 

```@repl quick-start
rydberg_populations = map(1:nsites) do i
    real(expect(put(nsites, i=>Op.n), prob.reg))
end
```
`prob.reg` is the register storing the final state after the time evolution.


## Looking for Help?

- check the slack channel [#julia](https://quera-workspace.slack.com/archives/C011C12GXRD)
- if not urgent, ask questions in [discussions](https://github.com/Happy-Diode/Bloqade.jl/discussions)

## Have Suggestions or Interested in Contribution?

- [file an issue](https://github.com/Happy-Diode/Bloqade.jl/issues/new) to report a bug or request a feature
