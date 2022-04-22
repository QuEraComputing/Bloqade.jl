# GPU Acceleration

Bloqade supports CUDA acceleration. To use
CUDA acceleration, you will need a NVIDIA graphics processing unit (GPU).

## Installation

To use CUDA accelerators, you need to install the CUDA package

```julia
pkg> add CUDA
```

This will automatically download all the needed dependencies of
the CUDA toolkit.

## Using CUDA

Converting your CPU-based simulation to CUDA-based simulation
is extremely simple: just use the `cu` function from `CUDA`
on the register object, which will convert the CPU-based
register to a CUDA-based register, e.g.

```julia
using CUDA
reg = zero_state(5)
dreg = cu(reg) # device register
```

For emulation, you can call `cu` on your emulation object
to convert everything (emulation intermediate memory, etc.)
into the GPU memory, e.g.

```julia
cu(KrylovEvolution(reg, clocks, h))
```

Other codes in Bloqade should work with CUDA
automatically.
