```@meta
CurrentModule = Bloqade
```

# Waveforms

Waveforms are essential ingredients for Rydberg quantum simulations. By controlling the waveforms of ``\Omega`` and ``\Delta``, one can prepare the ground states of certain target Hamiltonians and study their non-equalibrium dynamics. With Bloqade, we support several built-in waveforms and allow the users to specify waveforms by inputing functions. We also support different operations of waveforms, such as smoothing, waveform, and composing, and more. 

The generated waveforms can be directly used to build time-dependent Hamiltonians, please see [Hamiltonians](@ref) section for more details. 

## Creating Waveforms

In Bloqade, the waveforms are defined as [`Waveform`](@ref) object,
which is a composition of a callable object and a real number `duration`.

```@docs
BloqadeWaveforms.Waveform
```

Bloqade gives users the flexibility to specify general waveforms by inputing functions. The following code constructs a sinusoidal waveform with time duration of ``4 \pi \mu s``

```@example waveform
using Bloqade
using BloqadePlots: draw, draw!
using PythonCall
plt = pyimport("matplotlib.pyplot")
waveform = Waveform(t->2.2sin(t), duration=4π);
draw(waveform)
```
where `BloqadePlots` is a plotting package for objects from `Bloqade`,
that you need to use explicitly. And in our documentation we use the
python package [`matplotlib`](https://matplotlib.org) for plotting.

Bloqade supports built-in waveforms for convenience (see References below). 
For example, the codes below create different waveform shapes with single lines:

```@example waveform
waveform = piecewise_linear(clocks=[0.0, 0.2, 0.5, 0.8, 1.0], values=[0.0, 1.5, 3.1, 3.1, 0.0]); 
draw(waveform)
```

```@example waveform
waveform = piecewise_constant(clocks=[0.0, 0.2, 0.5, 0.7], values=[0.0, 1.5, 3.1]);
draw(waveform)
```

```@example waveform
waveform = linear_ramp(duration=0.5, start_value=0.0, stop_value=1.0);
draw(waveform)
```

```@example waveform
waveform =  constant(duration=0.5, value=2.1);
draw(waveform)
```

```@example waveform
waveform = sinusoidal(duration=4π, amplitude=2.2); 
draw(waveform)
```

In certain cases, users may have their own waveforms specified by a vector of clocks and a vector of signal strengths. To build a waveform from the two vectors, we can directly use the functions `piecewise_linear` or `piecewise_constant`, corresponding to different interpolations. 

```@example waveform
clocks = collect(0:1e-1:2);
values = rand(length(clocks));
wf1 = piecewise_linear(;clocks, values); 
wf2 = piecewise_constant(;clocks, values); 

fig, (ax1, ax2) = plt.subplots(figsize=(12, 4), ncols=2)
draw!(ax1, wf1)
draw!(ax2, wf2)
fig
```

For more advanced interpolation options, please see the [JuliaMath/Interpolations](http://juliamath.github.io/Interpolations.jl/latest/) package.

## Operations of Waveforms

Bloqade also supports several operations of the waveforms. 
Waveforms can be sliced using the duration syntax `start..stop`, e.g

```@example waveform
wf = sinusoidal(duration=2.2);
wf[1.1..1.5];
draw(wf)
```

Waveforms can be composed together via `append`

```@example waveform
wf1 = Waveform(sin, duration=2.2);
wf2 = linear_ramp(;start_value=0.0, stop_value=1.1, duration=0.5);
waveform = append(wf1, wf2); 
draw(waveform)
```

where the waveform `w2` is appended at the end of `w1`. 

Sharp waveforms may result in bad performance in practice (e.g. for adibatic preparing a ground state of target Hamiltnonian),
it is sometimes preferred to smoothen the waveform using
the moving average methods, one can use the [`smooth`](@ref)
function to create a smooth-ed wavefrom from a piecewise linear
waveform.

```@example waveform
wf = piecewise_linear(clocks=[0.0, 2.0, 3.0, 4.0], values=[0.0, 3.0, 1.1, 2.2]);
swf = smooth(wf);
draw(swf)
```

## Waveform Arithmetics

Bloqade also supports several arithmetics of waveforms. If two waveforms have the same duration, we can directly add up or subtract the strength of two waveforms, simply by using `+` or `-`. 

```@example waveform
wf1 = linear_ramp(;duration=2.2, start_value=0.0, stop_value=1.0);
wf2 = Waveform(sin, duration=2.2);
wf3 = wf1 + wf2; 
wf4 = wf1 - wf2;

fig, (ax1, ax2) = plt.subplots(figsize=(12, 4), ncols=2)
draw!(ax1, wf3)
draw!(ax2, wf4)
fig

```

If we want to increase the strength of a waveform by some times, we can directly use `*`

```@example waveform
wf = linear_ramp(;duration=2.2, start_value=0.0, stop_value=1.0);
wf_t = 3 * wf;

fig, (ax1, ax2) = plt.subplots(figsize=(12, 4), ncols=2)
draw!(ax1, wf)
draw!(ax2, wf_t)
fig

```

Such operation could also be broadcasted by using `.*`
```@example waveform
wf2, wf3 = [2.0, 3.0] .* wf1; 

fig, (ax1, ax2) = plt.subplots(figsize=(12, 4), ncols=2)
draw!(ax1, wf2)
draw!(ax2, wf3)
fig
```


## References


```@docs
piecewise_linear
piecewise_constant
linear_ramp
constant
sinusoidal
smooth
smooth(kernel, Xi::Vector, Yi::Vector, kernel_radius::Real)
```