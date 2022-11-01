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
ϕ = constant(;duration=T, value=0)
Ω = constant(;duration=T, value=0)

h = rydberg_h(atoms;Δ=Δi,ϕ=ϕ,Ω=Ω)
h_hardware, info = hardware_transform(h)
h = to_json(h_hardware)


open("lib/BloqadeSchema/schema_examples/aws_science_tasks/lattice_programmability/lattice3.json","w") do f
    JSON.print(f, JSON.parse(h))
end
