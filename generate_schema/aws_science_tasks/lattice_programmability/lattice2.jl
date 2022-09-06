using BloqadeSchema
using Bloqade
using JSON
# sites are spaced less than the minimum distance apart
atoms = AtomList([(0.0, 0.0), (0.0, 0.8), (0.0, 3.4), (0.0, 5.0)])

T=1
Δ = Waveform(t->cos(2π*t/T),T)
Δi=Δ
H = rydberg_h(atoms;Δ=Δi)
h = BloqadeSchema.to_json(H,waveform_tolerance=1e-1)



open("generate_schema/aws_science_tasks/lattice_programmability/lattice2.json","w") do f
    JSON.print(f, h)
end