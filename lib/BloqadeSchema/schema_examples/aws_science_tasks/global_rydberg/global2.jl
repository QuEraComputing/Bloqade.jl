using BloqadeSchema
using Bloqade
using JSON


# No detuning, and using a rabi_frequency pulse strength for a longer time, and see how much coherence is lost due to longer evolution time.

atoms = AtomList([(0.0, 0.0)])

T=2
# epsilon =0.01
# Ω = piecewise_constant(clocks=[0.0, epsilon, T+ epsilon , T+2*epsilon], values= 2π*[0.0, 3, 0])
Ω = constant(;duration=T, value=6π)
ϕ = constant(;duration=T, value=0)
Δ = constant(;duration=T, value=0)

h = rydberg_h(atoms;Ω=Ω,ϕ=ϕ,Δ=Δ)
h_hardware,info = hardware_transform(h)
h = to_json(h_hardware)


open("lib/BloqadeSchema/schema_examples/aws_science_tasks/global_rydberg/global2.json","w") do f
    JSON.print(f, JSON.parse(h))
end