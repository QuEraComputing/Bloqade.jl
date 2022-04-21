# CUDA Acceleration

The emulator supports CUDA acceleration, to use
CUDA acceleration, you will need NVIDIA graphic card.

## Installation

To use CUDA accelerators, you need to install CUDA package

```julia
pkg> add CUDA
```

This will automatically download all needed dependencies of
CUDA toolkit.

## Using CUDA

Converting your CPU-based simulation to CUDA-based simulation
is simple, just use the `cu` function from `CUDA`
on the register object, which will convert the CPU-based
register to a CUDA-based register, e.g

```julia
using CUDA
reg = zero_state(5)
dreg = cu(reg) # device register
```

For emulation, you can call `cu` on your emulation object
to convert everything (emulation intermediate memory etc.)
into GPU memory, e.g

```julia
cu(KrylovEvolution(reg, clocks, h))
```

other code in the emulator should adapt to CUDA
automatically.
