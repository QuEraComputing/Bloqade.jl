using BloqadeSchema
using Bloqade
using JSON


# 4-10 atoms positioned at the edge of the addressable area

atoms = AtomList([(0.0, 0.0), (0, 100), (100,99), (100, 0), (6, 7), (9, 6)])

T=1
Δ = Waveform(t->cos(2π*t/T),T)
Δi=Δ
Ω = Waveform(t->sin(2π*t/T)^2,T)
H = rydberg_h(atoms;Ω=Ω,Δ=Δi)
h = BloqadeSchema.to_json(H,waveform_tolerance=1e-1)


open("generate_schema/aws_science_tasks/lattice_programmability/lattice1.json","w") do f
    JSON.print(f, h)
end