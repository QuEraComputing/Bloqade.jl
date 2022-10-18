using BloqadeSchema
using Bloqade
using JSON

atom_coordinates = [
    [(i, 2) for i in range(5,9)],
    [(i, 3) for i in range(4,10)],
    [(3, 4), (4, 4), (5, 4), (9, 4), (10, 4), (11, 4)],
    [(3, 5), (4, 5), (10, 5), (11, 5)],
    [(3, 6), (4, 6), (10, 6), (11, 6)],
    [(9, 7), (10, 7), (11, 7)],
    [(8, 8), (9, 8), (10, 8)], 
    [(7, 9), (8, 9), (9, 9)], 
    [(7,10), (8, 10)],
    [(7,11), (8, 11)],
    [(7,13), (8, 13)],
    [(7,14), (8, 14)]
]

flattened_coordinates = vcat(atom_coordinates...)

# scale to satisfy minimum distance constraints (4.0 μm)
min_distance = 4.0
atoms = AtomList(map(c -> c .* min_distance, flattened_coordinates))

T=1/8

Ω = constant(;duration=T,value=8π)
ϕ = constant(;duration=T,value=0)
Δ = constant(;duration=T,value=0)

h = rydberg_h(atoms;Ω=Ω,ϕ=ϕ,Δ=Δ)
h_hardware,info = hardware_transform(h)
h = to_json(h_hardware)


open("lib/BloqadeSchema/schema_examples/question_mark.json","w") do f
    JSON.print(f, JSON.parse(h))
end
