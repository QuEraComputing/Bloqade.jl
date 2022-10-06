using BloqadeSchema
using Bloqade
using JSON


# 4-10 atoms positioned at the edge of the addressable area

atoms = AtomList([(0.0, 0.0), (0, 100), (100,99), (100, 0), (6, 7), (9, 6)])

T=1
Δ = Waveform(t->cos(2π*t/T),T)
Δi=Δ
Ω = Waveform(t->sin(2π*t/T)^2,T)
ϕ = constant(;duration=T,value=0)

h = rydberg_h(atoms;Ω=Ω,Δ=Δi,ϕ=ϕ)
h_hardware, info = hardware_transform(h)
h = to_json(h_hardware)

open("lib/BloqadeSchema/schema_examples/aws_science_tasks/lattice_programmability/lattice1.json","w") do f
    JSON.print(f, JSON.parse(h))
end