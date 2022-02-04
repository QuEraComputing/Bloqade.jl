# Quick Start

## Create a Hamiltonian

## Define Atom Positions

## Run emulation

## Create a Lattice

Creating a lattice is very simple in EaRyd, e.g we can create a square lattice as following

```@repl quick-start
using EaRyd
generate_sites(SquareLattice(), 3, 3)
```

this generates the atom positions on a ``3\times 3`` square lattice using the [`generate_sites`](@ref)
function.

we support the following built-in lattice: [`SquareLattice`](@ref), [`KagomeLattice`](@ref), [`HoneycombLattice`](@ref), and more. Please refer to [Lattices](@ref) for more detailed guide of lattice related operation.

!!! tips
    The lattice in EaRyd is actually defined as general Bravis lattice. You can create
    your own by defining the corresponding lattice vector. Check [Bravais Lattice](@ref bravais-lattice)

## Create a Waveform


EaRyd gives users the flexibility to specify general waveform by inputing functions. The following code constracting a sinusoidal waveform with time duration of ``4 \pi``

```julia
waveform = Waveform(t->2.2sin(t), duration=4π)
```

We also support several built-in time-dependent waveforms, including [`piecewise_linear`](@ref), [`piecewise_constant`](@ref), [`linear_ramp`](@ref), [`constant`](@ref), [`sinusoidal`](@ref). For example, we can create a piecewise linear waveform simply by one line below 

```@repl quick-start
waveform = piecewise_linear(clocks=[0.0, 0.2, 0.5, 0.8, 1.0], values=[0.0, 1.5, 3.1, 3.1, 0.0])
```

Please refer to [Waveform](@ref) for more detailed guide of waveform related operation.


## Create a Register