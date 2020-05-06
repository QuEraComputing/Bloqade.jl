using BenchmarkTools
using RydbergEmulator
using LightGraphs
using Random

Random.seed!(42)
n = 32; depth = 5
g = SimpleGraph(n)
for i in 1:n, j in 1:n
    if rand() < 0.2
        add_edge!(g, i, j)
    end
end

subspace_v = subspace(g)
hs = SimpleRydberg.(rand(depth))
ts = rand(depth)
r = RydbergEmulator.zero_state(n, subspace_v)
qaoa = QAOA{n}(subspace_v, hs, ts)

@benchmark r |> qaoa

@profiler r |> qaoa
