```@meta
CurrentModule = EaRydWaveforms
```

# Waveform

Waveforms are signals used in pulse programming or as hamiltonian parameters. The waveforms are defined as [`Waveform`](@ref) object,
which is a composition of a callable object and a real number `duration`.

## Creating Waveforms

EaRyd gives users the flexibility to specify general waveform by inputing functions. The following code constracting a sinusoidal waveform with time duration of ``4 \pi``

```@repl creating-waveform
using EaRydWaveforms # hide
waveform = Waveform(t->2.2sin(t), duration=4π)
```

We also support several built-in time-dependent waveforms, including [`piecewise_linear`](@ref), [`piecewise_constant`](@ref), [`linear_ramp`](@ref), [`constant`](@ref), [`sinusoidal`](@ref). For example, we can create a piecewise linear waveform simply by one line below 

```@repl creating-waveform
waveform = piecewise_linear(clocks=[0.0, 0.2, 0.5, 0.8, 1.0], values=[0.0, 1.5, 3.1, 3.1, 0.0])
```

## Slicing Waveforms

Waveforms can be sliced using the duration syntax `start..stop`, e.g

```@repl slicing
using EaRydWaveform # hide
wf = sinusoidal(duration=2.2)
wf[1.1..1.5]
```

## Composing Waveforms

Waveforms can be composed together via `append`

```@docs
EaRydWaveforms.append
```

## Waveform Smoothing

Sharp waveforms may result in bad performance in practice,
it is sometimes preferred to smoothen the waveform using
the moving average methods, one can use the [`smooth`](@ref)
function to create a smooth-ed wavefrom from a piecewise linear
waveform.
