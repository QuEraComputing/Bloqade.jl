using Test
using RydbergEmulator
using LightGraphs

g = SimpleGraph(4)
add_edge!(g, 1, 2)
add_edge!(g, 1, 3)
add_edge!(g, 1, 4)
add_edge!(g, 2, 3)
add_edge!(g, 2, 4)
add_edge!(g, 3, 4)

Ω = Float64[1, 1, 1, 1]
Δ = Float64[1, 1, 1, 1]
ϕ = Float64[π, π, π, π]

to_matrix(g, Ω, ϕ, Δ)
