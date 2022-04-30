# Working with Units

Unit conversion can be both tedious and prone to errors when working with simulation
of real quantum hardware. 
In Bloqade, the physical variables can support units explicitly from the [Unitful](https://github.com/PainterQubits/Unitful.jl) package.
Explicit units will be converted to our
default units set automatically.  In other words,  you can use a different unit set to create the Hamiltonian other than the default units, and not worry about conversion yourself. 
For detailed information about 
the default unit of the Hamiltonian, please refer to the [Bloqade](@ref) page.

```@repl units
using Bloqade
using Unitful: kHz, µm
rydberg_h([(1, ), (2, )], C = 2π * 109.2kHz * µm^6)
```
In the above example, we have assigned the parameter `C_6` with the unit `kHz * μm^6`, which is automatically
converted to the default unit `MHz * µm^6` in Bloqade. 

Other than the Hamiltonian, we can also specify units on waveforms, e.g.: 

```@example units
using Unitful: rad, ms

wf = piecewise_linear(clocks=[0.0ms, 0.1ms, 0.2ms], values= [0.1, 1.1, 2.1] .* (rad/ms))
Bloqade.plot(wf)
```
From this plot, we can see that the units for `clocks` and `values` have been converted to `μs` and `MHz` respectively. 
