# [GPU Acceleration] (@id cuda)

!!! danger "CUDA.jl Version 3 Required"

    Bloqade currently does not work with the latest version (4.x) of the [CUDA](https://github.com/JuliaGPU) package. As a workaround you will need to install version 3. The installation instructions below should already address this.

Bloqade supports CUDA acceleration. To use
CUDA acceleration, you will need a NVIDIA graphics processing unit (GPU).

## Installation

To use CUDA accelerators, you need to install the [CUDA](https://github.com/JuliaGPU/CUDA.jl) package and the [Adapt](https://github.com/JuliaGPU/Adapt.jl) package:

```julia
pkg> add CUDA@3 Adapt
```

This will automatically download all the needed dependencies of
the CUDA toolkit. This functionality requires CUDA toolkit
11.6+ and CUSPARSE 11.4+. You can check your version via

```julia
julia> CUDA.version()
```

## Using CUDA

Converting your CPU-based simulation to CUDA-based simulation
is extremely simple: just use the `adapt` function from `Adapt`
on the register object, which will convert the CPU-based
register to a CUDA-based register, e.g.:

```julia
using CUDA, Adapt
reg = zero_state(5)
dreg = adapt(CuArray, reg) # device register
```

For emulation, you can call `adapt` on your emulation object
to convert everything (emulation intermediate memory, etc.)
into the GPU memory, e.g.:

```julia
adapt(CuArray, KrylovEvolution(reg, clocks, h))
```

To perform operations after emulation (such as calculating the rydberg density) 
you will need to move the object out of GPU memory. This can be done so via `adapt` again:

```julia
# get the problem onto the GPU and perform emulation there:
gpu_prob = adapt(CuArray, SchrodingerProblem(register, time_span, hamiltonian))
integrator = init(gpu_prob, Vern9())
emulate!(integrator)

# get the problem off the GPU for subsequent processing/analysis:
prob = adapt(Array, gpu_prob)
```

Other code written with Bloqade should work with CUDA
automatically.
