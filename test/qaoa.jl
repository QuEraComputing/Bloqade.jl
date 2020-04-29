using Test
using Yao
using RydbergEmulator
using LightGraphs
using ExponentialUtilities
using LinearAlgebra
using RydbergEmulator: subspace

function naive_qaoa(st, g, hs, ts)
    for (h, t) in zip(hs, ts)
        st = expv(-im * t, to_matrix(g, 1.0, h.ϕ), st)
    end
    return st
end

@testset "standard qaoa" begin
    g = SimpleGraph(5)
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)
    add_edge!(g, 2, 4)
    add_edge!(g, 2, 5)
    add_edge!(g, 3, 4)
    add_edge!(g, 4, 5)


    hs = SimpleRydberg.(rand(10))
    ts = rand(10)
    subspace_v = [0, 1, 2, 4, 5, 8, 9, 16, 17, 20, 21]
    r = RydbergEmulator.zero_state(5, subspace_v)
    qaoa = QAOA{5}(subspace_v, hs, ts)
    target = naive_qaoa(copy(vec(r.state)), g, hs, ts)
    r |> qaoa
    @test  vec(r.state) ≈ target
    @test isnormalized(r)

    new_hs = rand(10)
    new_ts = rand(10)

    update_ansatz!(qaoa, new_hs, new_ts)
    @test map(x->x.ϕ, qaoa.hamiltonians) ≈ new_hs
    @test qaoa.ts ≈ ts
end
