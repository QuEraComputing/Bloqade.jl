using Test
using RydbergEmulator
using LightGraphs

g = SimpleGraph(5)
add_edge!(g, 1, 2)
add_edge!(g, 2, 3)
add_edge!(g, 2, 4)
add_edge!(g, 2, 5)
add_edge!(g, 3, 4)
add_edge!(g, 4, 5)

hs = simple_rydberg.(rand(10))
ts = rand(10)
subspace_v = [0, 1, 2, 4, 5, 8, 9, 16, 17, 20, 21]
s = Subspace(subspace_v)
r = RydbergEmulator.zero_state(5, subspace_v)
qaoa = QAOA{5}(ts, hs; subspace=s)

r |> qaoa
