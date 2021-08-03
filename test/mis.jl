using Test
using RydbergEmulator
using BitBasis

# generate random atom positions
atoms = RydAtom.([(0.0, 1.0), (1.0, 0.), (2.0, 0.0),
(1.0, 1.0), (1.0, 2.0), (2.0, 2.0)])
graph = unit_disk_graph(atoms, 1.5)
config = [1, 1, 1, 0, 1, 1]

@test RydbergEmulator.num_mis_violation(config, graph, 1) == 2
@test RydbergEmulator.num_mis_violation(config, graph, 2) == 2
@test RydbergEmulator.num_mis_violation(config, graph, 3) == 1
@test RydbergEmulator.num_mis_violation(config, graph, 4) == 0
@test RydbergEmulator.num_mis_violation(config, graph, 5) == 2
@test RydbergEmulator.num_mis_violation(config, graph, 6) == 1

@test !is_independent_set(graph, config)
to_independent_set!(graph, config)
@test is_independent_set(graph, config)


config = bitarray(36, length(atoms))
RydbergEmulator.to_independent_set!(graph, config)

space = blockade_subspace(graph)
raw_state = zeros(ComplexF64, length(space))
raw_state[space.map[packbits(config)]] = 1.0
r = RydbergReg{length(atoms)}(raw_state, space)
@test reduced_mean_rydberg(r, graph) == mean_rydberg(r)

# TODO: add an violation test
