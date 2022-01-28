# Quick Start

In this quick start, we will use the Rydberg Emulator to build a interacting Rydberg Hamiltonian, described by

$\frac{\mathcal{H}}{\hbar} = \sum_i \frac{\Omega_i}{2} \sigma_x^i - \sum_i \Delta_i n_i + \sum_{i < j} V_{ij} n_i n_j$

where $\Delta_i$ are the detunings of the driving lasers from the Rydberg state, $\sigma_x^i = |g_i\rangle \langle r_i| + |r_i\rangle \langle g_i|$ describes the coupling between the ground state $|g_i\rangle$ and the Rydberg state $|r_i\rangle$ of an atom at position $i$, driven at Rabi frequencey $\Omega_i$, $n_i = |r_i\rangle \langle r_i|$, $V_{ij} = C/r_{ij}^6$ is the interactions between atom $i$ and $j$, $C$ is the interacting constant that depends on particular Rydberg atoms, and $\hbar$ is the reduced Plank's constant. 


We will also show how to emulate quantum many-body dynamics governed by such a Hamiltonian.  

We start by loading the Emulator Module

```@repl quick-start
using EaRyd
```
## Define Atom Positions

As one can see from the Rydberg Hamiltonian, the interacting strengths between Rydberg atoms depend on their positions. To allow the users easily generate atom positions, EaRyd provides several built-in lattice structures. For instance, we can use the following codes to generate a chain of 10 atoms in 1D 

```@repl quick-start
nsites = 10
atoms = generate_sites(ChainLattice(), nsites, scale=5.72)
```
Here by specifying the scale, we setting the distance between nearest neighbor atoms to be 5.72 $\mu m$. Note that the default unit of length in EaRyd is $\mu m$. Correspondingly, the default value for interacting constant $C$ is $2π * 858386$ MHz $\times \mu m^6$ for Rydberg atom $^{87}Rb$ and $70 s$ Rydberg state. 

## Create a Hamiltonian

To create the full Hamiltonian, we will also need to specify the parameters $\Omega$ and $\Delta$. EaRyd supports both time independent and time dependent waveforms for $\Omega$ and $\Delta$. Several built-in waveforms and waveform functions are also supported. Here we will focus on a simple example where both $\Omega$ and $\Delta$ are constant with time, 
```@repl quick-start
Ω0=4π
Δ0=0
```
The default energy unit (also for $\Omega$ and $\Delta$) is MHz in EaRyd. Now that all the parameters specified, we can create a interacting Rydberg Hamiltonian, 

```@repl quick-start
h = rydberg_h(atoms;C = 2π * 858386, Ω=Ω0, Δ=Δ0)
```
 The above Hamiltonian is a global Hamiltonian without site-dependence because $\Omega$ and $\Delta$ are both single numbers. By making them as vectors, one can also create a site-dependent Hamiltonian. 


## Create a Register

We use the following code to create an initial state with all the atoms in the ground state 

```@repl quick-start
init= zero_state(10)
```

## Run emulation

We first create the problem where we want to emulate the system 
starting from the intial state and evolving under the defined Hamiltoninan. Here we choose the ODE solver.  

```@repl quick-start
prob = ODEEvolution(init, 1.6, h)
```
Here 1.6 is the total evolution time and its default unit is $\mu s$. Then we use the following code to emulate the problem. 

```@repl quick-start
emulate!(prob)
```

After the time-evolution, we are able to measure the expectation value of Rydberg density operator for each site. 

```@repl quick-start
densities = map(1:nsites) do i
    real(expect(put(nsites, i=>Op.n), prob.reg))
end
```
Here the prob.reg is a register storing the final state after time-evolution, and put(nsites, i=>Op.n) specifies the observable to be Rydberg density operator on site i. 
