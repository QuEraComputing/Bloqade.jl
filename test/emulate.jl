using Test
using Yao
using SparseArrays
using RydbergEmulator
using LightGraphs
using LinearAlgebra

if !isdefined(@__MODULE__, :test_graph)
    include("utils.jl")
end

function simple_evolve!(r::AbstractRegister, ts, hs, s::Subspace)
    st = vec(r.state)
    for (h, t) in zip(hs, ts)
        st = exp(-im * t * Matrix(SparseMatrixCSC(h, s))) * st
        normalize!(st)
    end
    r.state .= st
    return r
end

function simple_evolve!(r::AbstractRegister, ts, hs)
    st = vec(r.state)
    for (h, t) in zip(hs, ts)
        st = exp(-im * t * Matrix(SparseMatrixCSC(h))) * st
    end
    r.state .= st
    return r
end

@testset "subspace qaoa" begin
    hs = simple_rydberg.(5, rand(5))
    ts = rand(5)
    s = Subspace(test_subspace_v)
    r = RydbergEmulator.zero_state(5, s)

    cache = DiscreteEmulationCache(Complex{eltype(ts)}, first(hs), s)
    r1 = emulate!(copy(r), ts, hs; cache=cache)
    r2 = simple_evolve!(copy(r), ts, hs, s)

    @test r1 ≈ r2
end

@testset "fullspace" begin
    hs = simple_rydberg.(4, rand(4))
    ts = rand(4)
    r = Yao.zero_state(4)
    cache = DiscreteEmulationCache(Complex{eltype(ts)}, first(hs))
    r1 = emulate!(copy(r), ts, hs; cache=cache)
    r2 = simple_evolve!(copy(r), ts, hs)
    @test r1 ≈ r2
end

@testset "emulate" begin
    hs = simple_rydberg.(5, rand(3))
    ts = rand(3)
    @test emulate(ts, hs) ≈ emulate!(Yao.zero_state(5), ts, hs)
    @test emulate(test_subspace, ts, hs) ≈ emulate!(RydbergEmulator.zero_state(5, test_subspace), ts, hs)
end

@testset "trotterize" begin
    atoms = square_lattice(10, 0.8)
    space = blockade_subspace(atoms, 1.5)
    h = rydberg_h(atoms, 1.0, 2.0, sin)
    prob = DiscreteEvolution(zero_state(length(atoms), space), 0.5, h)
    display(prob)

    @test prob.t_or_ts ≈ 0:1e-3:0.5
    emulate!(prob)
    ts = [1e-3 for _ in 0:1e-3:0.5]
    hs = [h(t) for t in 0:1e-3:0.5]
    target = simple_evolve!(zero_state(length(atoms), space), ts, hs, space)
    @test prob.reg ≈ target

    prob = DiscreteEvolution(zero_state(length(atoms), space), 0.001, h; progress=true)
    @test_logs (:info, "emulating") (:info, "emulating") emulate!(prob)
    @test prob.options.progress == true
end
