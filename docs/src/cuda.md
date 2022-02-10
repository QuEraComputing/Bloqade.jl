# CUDA Acceleration

The emulator supports CUDA acceleration, to use
CUDA acceleration, you will need NVIDIA graphic card.

Converting your CPU-based simulation to CUDA-based simulation
is simple, just use the `cu` function from `CUDA`
on the register object, which will convert the CPU-based
register to a CUDA-based register, e.g

```julia
reg = zero_state(5)
dreg = cu(reg) # device register
```

other code in the emulator should adapt to CUDA
automatically.
