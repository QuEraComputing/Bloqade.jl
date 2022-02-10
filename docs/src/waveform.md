```@meta
CurrentModule = EaRydWaveforms
```

# Waveform

Waveforms are signals used in pulse programming or as hamiltonian parameters. The waveforms are defined as [`Waveform`](@ref) object,
which is a composition of a callable object and a real number `duration`.

```@docs
EaRydWaveforms.Waveform
```

## Creating Waveforms

The following builtin waveforms are provided:

```@docs
piecewise_linear
piecewise_constant
linear_ramp
constant
sinusoidal
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
the moving average methods, we provide `smooth` function
to create a smooth-ed wavefrom from a piecewise linear
waveform.

```@docs
EaRydWaveforms.smooth
```
