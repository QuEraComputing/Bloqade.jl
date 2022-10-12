using BloqadeSchema
using Bloqade
using JSON

# sites are spaced less than the minimum distance apart
atoms = AtomList([(0.0, 0.0), (0.0, 0.8), (0.0, 3.4), (0.0, 5.0)])

T=1
Δ = Waveform(t->cos(2π*t/T),T)
Δi=Δ
Ω = constant(;duration=T, value=0)
ϕ = constant(;duration=T, value=0)

h = rydberg_h(atoms;Δ=Δi,Ω=Ω,ϕ=ϕ)
h_hardware,info = hardware_transform(h)
h = to_json(h_hardware)


open("lib/BloqadeSchema/schema_examples/aws_science_tasks/lattice_programmability/lattice2.json","w") do f
    JSON.print(f, JSON.parse(h))
end