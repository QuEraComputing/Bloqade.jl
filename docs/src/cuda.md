# [GPU Acceleration] (@id cuda)

Bloqade supports CUDA acceleration. To use
CUDA acceleration, you will need a NVIDIA graphics processing unit (GPU).

## Installation

To use CUDA accelerators, you need to install the CUDA package and the Adapt package:

```julia
pkg> add CUDA Adapt
```

This will automatically download all the needed dependencies of
the CUDA toolkit. This functionality requires CUDA toolkit
11.6+ and CUSPARSE 11.4+. You can check your version via

```julia
julia> CUDA.version()
```

## Using CUDA

Converting your CPU-based simulation to CUDA-based simulation
is extremely simple: just use the `cu` function from `CUDA`
on the register object, which will convert the CPU-based
register to a CUDA-based register, e.g.:

```julia
using CUDA, Adapt
reg = zero_state(5)
dreg = adapt(CuArray, reg) # device register
```

For emulation, you can call `cu` on your emulation object
to convert everything (emulation intermediate memory, etc.)
into the GPU memory, e.g.:

```julia
adapt(CuArray, KrylovEvolution(reg, clocks, h))
```

Other codes in Bloqade should work with CUDA
automatically.
