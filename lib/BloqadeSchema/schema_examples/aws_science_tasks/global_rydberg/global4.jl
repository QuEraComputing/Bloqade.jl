
using BloqadeSchema
using Bloqade
using JSON

# Switching rapidly between 0 and non-zero detuning, while also switching rapidly, in sync, between maximal and zero rabi frequency.
# The result should be identical to applying a continuous rabi pulse over a timeframe when the detuning was zero.

atoms = AtomList([(0.0, 0.0)])

T1=1
T2=2.5
Ω = piecewise_constant(clocks=[0.0, T1, T2], values= 2π*[0.0, 4])

# maximal detuning 
Δ = piecewise_constant(clocks=[0.0, T1, T2], values= 2π*[20, 0])

ϕ = constant(;duration=T2, value=0)

h = rydberg_h(atoms;Ω=Ω, Δ=Δ, ϕ=ϕ)

# hardware_transform will add slopes between values in waveform when necessary
# with the maximum slew rate supported by hardware
h_hardware,info  = hardware_transform(h)
h = to_json(h_hardware)

open("lib/BloqadeSchema/schema_examples/aws_science_tasks/global_rydberg/global4.json","w") do f
    JSON.print(f, JSON.parse(h))
end

