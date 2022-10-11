using BloqadeSchema
using Bloqade
using JSON


# No detuning, and constant rabi_frequency pulse applied at max strength for a short period of time (to achieve a quick pi-pulse)


atoms = AtomList([(0.0, 0.0)])

T=1/8

Ω = constant(;duration=T,value=8π)
ϕ = constant(;duration=T,value=0)
Δ = constant(;duration=T,value=0)

h = rydberg_h(atoms;Ω=Ω,ϕ=ϕ,Δ=Δ)
h_hardware,info = hardware_transform(h)
h = to_json(h_hardware)


open("lib/BloqadeSchema/schema_examples/aws_science_tasks/global_rydberg/global1.json","w") do f
    JSON.print(f, JSON.parse(h))
end

