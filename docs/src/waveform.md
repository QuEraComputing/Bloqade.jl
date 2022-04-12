```@meta
CurrentModule = BloqadeWaveforms
```

# Waveforms

Waveforms are essential ingredients for Rydberg quantum simulations. By controlling the waveforms of ``\Omega`` and ``\Delta``, one can prepare ground states of certain target Hamiltonian and study non-equalibrium dynamics of time-dependent Hamiltonians. With Bloqade, we support several built-in waveforms and allow the users to specify waveforms by inputing functions. We also support different operations of waveforms, such as smoothing, slicing, and composing, et al. 

The generated waveforms can be directly used to build time-dependent Hamiltonians, see [`Hamiltonian`](@ref). 

## Creating Waveforms

In Bloqade, the waveforms are defined as [`Waveform`](@ref) object,
which is a composition of a callable object and a real number `duration`.

Bloqade gives users the flexibility to specify general waveforms by inputing functions. The following code constracting a sinusoidal waveform with time duration of ``4 \pi``

```@repl creating-waveform
using Bloqade
using BloqadePlots
waveform = Waveform(t->2.2sin(t), duration=4π);
draw(waveform)
```

We also support several built-in time-dependent waveforms, including [`piecewise_linear`](@ref), [`piecewise_constant`](@ref), [`linear_ramp`](@ref), [`constant`](@ref), [`sinusoidal`](@ref). For example, the codes below create different waveform shapes simply by single lines 

```@repl creating-waveform
waveform = piecewise_linear(clocks=[0.0, 0.2, 0.5, 0.8, 1.0], values=[0.0, 1.5, 3.1, 3.1, 0.0]); 
draw(waveform)
```

```@repl creating-waveform
waveform = piecewise_constant(clocks=[0.0, 0.2, 0.5], values=[0.0, 1.5, 3.1]); 
draw(waveform)
```

```@repl creating-waveform
waveform = linear_ramp(duration=0.5, start_value=0.0, stop_value=1.0);
draw(waveform)
```

```@repl creating-waveform
waveform =  constant(duration=0.5, value=2.1);
draw(waveform)
```

```@repl creating-waveform
waveform =  sinusoidal(duration=4π, amplitude=2.2); 
draw(waveform)
```

In certain cases, users may have their own waveforms specified by a vector of clocks and a vector of signal strengths. To build a waveform from the two vectors, we can directly use the functions `piecewise_linear` or `piecewise_constant`, corresponding to different interpolations. 

```@repl user_input
a =  [i for i in 1:10]
b = rand(10)
wf1 = piecewise_linear(clocks=a; values=b)
wf2 = piecewise_constant(clocks=a; values=b)
```

## Operations of Waveforms

Bloqade also supports several operations of the waveforms. 
Waveforms can be sliced using the duration syntax `start..stop`, e.g

```@repl slicing
using BloqadeWaveform # hide
wf = sinusoidal(duration=2.2);
wf[1.1..1.5];
draw(wf)
```

Waveforms can be composed together via `append`
```@repl append
wf1 = Waveform(sin, duration=2.2)
wf2 = linear_ramp(;start_value=0.0, stop_value=1.1, duration=0.5)
waveform = append(wf1, wf2)
```
where the waveform `w2` is appended at the end of `w1`. 


Sharp waveforms may result in bad performance in practice (e.g. for adibatic preparing a ground state of target Hamiltnonian),
it is sometimes preferred to smoothen the waveform using
the moving average methods, one can use the [`smooth`](@ref)
function to create a smooth-ed wavefrom from a piecewise linear
waveform.

```@repl append
wf = piecewise_linear(clocks=[0.0, 2.0, 3.0, 4.0], values=[0.0, 3.0, 1.1, 2.2])
swf = smooth(wf)
```

## Waveform Linear Algebra

Bloqade also supports several linear algebra operations of waveforms. If two waveforms have the same duration, we can directly add up or subtract the strength of two waveforms, simply by using `+` or `-`. 

```@repl add
wf1 = linear_ramp(;duration=2.2, start_value=0.0, stop_value=1.0);
wf2 = Waveform(sin, duration=2.2);
wf3 = wf1 + wf2
wf4 = wf1 - wf2
```

If we want to increase the strength of a waveform by some times, we can directly use `*`

```@repl times
wf = linear_ramp(;duration=2.2, start_value=0.0, stop_value=1.0);
wf_t = 3 * wf
```

Such operation could also be broadcasted by using `.*`
```@repl times
wf2, wf3 = [2.0, 3.0] .* wf1
```

