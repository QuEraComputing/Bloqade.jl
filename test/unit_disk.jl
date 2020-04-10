using Test
using RydbergEmulator
using LightGraphs
atoms = AtomPosition(10, 2)

g = unit_disk_graph(atoms, 1.0)

for each in edges(g)
    @show each
end
