using BenchmarkTools
using RydbergEmulator
using LightGraphs
using Random

Random.seed!(42)
n = 10; depth = 5
g = SimpleGraph(n)
add_edge!(g, 1, 2)
add_edge!(g, 2, 3)
add_edge!(g, 2, 4)
add_edge!(g, 2, 5)
add_edge!(g, 3, 4)
add_edge!(g, 4, 5)
add_edge!(g, 5, 10)
add_edge!(g, 8, 5)
add_edge!(g, 7, 5)
add_edge!(g, 2, 10)

subspace_v = subspace(g)
hs = SimpleRydberg.(rand(depth))
ts = rand(depth)
r = RydbergEmulator.zero_state(n, subspace_v)
qaoa = QAOA{n}(subspace_v, hs, ts)

@benchmark r |> qaoa

@profiler r |> qaoa
