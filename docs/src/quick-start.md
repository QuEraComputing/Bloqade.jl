# Quick Start

In this quick start, we will use the Rydberg Emulator to build a interacting Rydberg Hamiltonian, described by

```math
\frac{\mathcal{H}}{\hbar} = \sum_i \frac{\Omega_i}{2} \sigma_x^i - \sum_i \Delta_i n_i + \sum_{i < j} V_{ij} n_i n_j
```

where ``\Delta_i`` are the detunings of the driving lasers from the Rydberg state, ``\sigma_x^i = |g_i\rangle \langle r_i| + |r_i\rangle \langle g_i|`` describes the coupling between the ground state ``|g_i\rangle`` and the Rydberg state ``|r_i\rangle`` of an atom at position ``i``, driven at Rabi frequencey ``\Omega_i``, ``n_i = |r_i\rangle \langle r_i|``, ``V_{ij} = C/r_{ij}^6`` is the interactions between atom ``i`` and ``j``, ``C`` is the interacting constant that depends on particular Rydberg atoms, and ``\hbar`` is the reduced Plank's constant. 


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
Here by specifying the scale, we setting the distance between nearest neighbor atoms to be 5.72 ``\mu m``. Note that the default unit of length in EaRyd is ``\mu m``. Correspondingly, the default value for interacting constant ``C`` is ``2π * 858386`` MHz ``\times \mu m^6`` for Rydberg atom ``^{87}Rb`` and ``70 s`` Rydberg state. 

## Create a Hamiltonian

To create the full Hamiltonian, we will also need to specify the parameters ``\Omega`` and ``\Delta``. EaRyd supports both time independent and time dependent waveforms for ``\Omega`` and ``\Delta``. Several built-in waveforms and waveform functions are also supported. Here we will focus on a simple example where both ``\Omega`` and ``\Delta`` are constant with time, 
```@repl quick-start
Ω0=4π
Δ0=0
```
The default energy unit (also for ``\Omega`` and ``\Delta``) is MHz in EaRyd. Now that all the parameters specified, we can create a interacting Rydberg Hamiltonian, 

```@repl quick-start
h = rydberg_h(atoms;C = 2π * 858386, Ω=Ω0, Δ=Δ0)
```
 The above Hamiltonian is a global Hamiltonian without site-dependence because ``\Omega`` and ``\Delta`` are both single numbers. By making them as vectors, one can also create a site-dependent Hamiltonian. 


## Create a Register

We use the following code to create an initial state with all the atoms in the ground state 

```@repl quick-start
init= zero_state(10)
```

## Run emulation

We are interested in calculating the real-time dynamics of Rydberg density for each site under the evolution of the defined Hamiltonian. The  total time step of the problem is choosen to be 120, and each time step being `` 0.01 \mu s``. Note that the default unit of time is ``\mu s`` in EaRyd. 

```@repl quick-start
iteration = 1:120
ts = [0.01 for _ in iteration];
hs = [h for _ in iteration];
clocks = cumsum(ts)
``` 

We can set up the problem by specifying the intial state, evolving time (vector), and the time-independent Hamiltoninan (vector) defined above. The  Krylov solver can be used for this example.  

```@repl quick-start
prob1 = KrylovEvolution(init, ts, hs)
```

The following for loop is used to calculate the expectation value of the Rydberg density for individual site and for each time step.  

```@repl quick-start 
density_site = zeros(nsites, length(iteration)); 
for info in prob1
    for i in 1:nsites
        density_site[i, info.step] = expect(put(nsites, i=>Op.n), info.reg)
    end
end
```
info.reg is a register storing the quantum state at each time step, and put(nsites, i=>Op.n) specifies the observable to be Rydberg density operator on site ``i``.


On the other hand, if we are only interested in calculating the expectation value at the final time step, we can directly emulate the problem first and then the measurement for the final state. 

```@repl quick-start
prob2 = ODEEvolution(init, 1.6, h)
emulate!(prob2)
```
Here we have choosen the ODE solver and the total evolution time is set to be 1.6 ``\mu s``
```@repl quick-start
densities = map(1:nsites) do i
    real(expect(put(nsites, i=>Op.n), prob2.reg))
end
```
prob.reg is the register storing the final state after the time-evolution. 