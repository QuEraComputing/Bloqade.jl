# Working with Units

Unit conversion can be tedious when working with simulation
of real quantum experiments. In Bloqade, the physical variables can support explicit units from [Unitful](https://github.com/PainterQubits/Unitful.jl).
When using explicit units, they will be converted to our
default unit set automatically.  In other words,  we can use a different unit set to create the Hamiltonian other than the default units. For detailed information about 
the detault unit of Hamiltonian, please refer to [Bloqade](@ref). 

```@repl units
using Bloqade
using Unitful: kHz, µm
rydberg_h([(1, ), (2, )], C=109.2kHz * µm^6)
```
In the above example, we have assigned the the parameter `C` with the unit `kHz * µm^6`, which is automatically
converted to default unit `MHz * µm^6` in Bloqade. 

Other than the Hamiltonian, we can also specify units on waveforms, e.g. 

```@repl units
using Unitful: rad, ms
wf = piecewise_linear(clocks=[0.0ms, 0.1ms, 0.2ms], values=[0.1, 1.1, 2.1] .* (rad/ms))
```
from the plot, we see that the units for `clocks` and `values` have been converted to `μs` and `MHz` respectively. 