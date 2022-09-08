using BloqadeSchema
using Bloqade
using JSON

# max number of atoms (256) arranged in various configurations

# square lattice 
atoms = generate_sites(SquareLattice(), 16, 16; scale = 4.5)

# Lieb lattice 
# atoms = generate_sites(LiebLattice(), 10, 8; scale = 4.5)

# TriangularLattice
# atoms = generate_sites(TriangularLattice(), 16, 16; scale = 3.2)



T=1
Δ = Waveform(t->cos(2π*t/T),T)
Δi=Δ
H = rydberg_h(atoms;Δ=Δi)
h = BloqadeSchema.to_json(H,waveform_tolerance=1e-1)



open("generate_schema/aws_science_tasks/lattice_programmability/lattice3.json","w") do f
    JSON.print(f, h)
end
