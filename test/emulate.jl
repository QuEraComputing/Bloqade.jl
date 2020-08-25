using Test
using Yao
using SparseArrays
using RydbergEmulator
using LightGraphs
using ExponentialUtilities
using LinearAlgebra

if !isdefined(@__MODULE__, :test_graph)
    include("utils.jl")
end

function naive_qaoa!(r::AbstractRegister, hs, ts, s::Subspace)
    st = vec(r.state)
    for (h, t) in zip(hs, ts)
        st = expv(-im * t, SparseMatrixCSC(h, s), st)
    end
    r.state .= st
    return r
end

function naive_qaoa!(r::AbstractRegister, hs, ts)
    st = vec(r.state)
    for (h, t) in zip(hs, ts)
        st = expv(-im * t, SparseMatrixCSC(h), st)
    end
    r.state .= st
    return r
end

@testset "subspace qaoa" begin
    hs = simple_rydberg.(5, rand(5))
    ts = rand(5)
    s = Subspace(test_subspace_v)
    r = RydbergEmulator.zero_state(5, s)

    cache = EmulatorCache(eltype(ts), first(hs), s)
    r1 = emulate!(copy(r), ts, hs; cache=cache)
    r2 = naive_qaoa!(copy(r), hs, ts, s)

    @test r1 ≈ r2
end

@testset "fullspace" begin
    hs = simple_rydberg.(4, rand(4))
    ts = rand(4)
    r = Yao.zero_state(4)
    cache = EmulatorCache(eltype(ts), first(hs), 4)
    r1 = emulate!(copy(r), ts, hs; cache=cache)
    r2 = naive_qaoa!(copy(r), hs, ts)
    @test r1 ≈ r2
end

@testset "contiguous time" begin
    h = XTerm(5, 1.0, sin)
    dt = 1e-5
    @testset "subspace" begin
        r1 = RydbergEmulator.zero_state(5, test_subspace)
        emulate!(r1, 0.2, h)

        r2 = RydbergEmulator.zero_state(5, test_subspace)
        emulate!(r2, map(_->dt, 0.0:dt:0.2), map(h, 0.0:dt:0.2))

        @test isapprox(r1.state, r2.state; atol=1e-4)
    end

    @testset "fullspace" begin
        r1 = Yao.zero_state(5)
        emulate!(r1, 0.2, h)
        r2 = Yao.zero_state(5)
        emulate!(r2, map(_->dt, 0.0:dt:0.2), map(h, 0.0:dt:0.2))
        @test isapprox(r1.state, r2.state; atol=1e-4)
    end
end
