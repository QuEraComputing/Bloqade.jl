using Test
using Yao
using SparseArrays
using RydbergEmulator
using LightGraphs
using ExponentialUtilities
using LinearAlgebra
using CUDA

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

    qaoa = QAOA(hs[1], s)
    r1 = copy(r) |> qaoa(ts, hs)
    r2 = naive_qaoa!(copy(r), hs, ts, s)

    @test r1 ≈ r2

    @testset "cuda" begin
        if CUDA.functional()
            dqaoa = cu(qaoa)
            @test isapprox((cu(r) |> dqaoa(ts, hs) |> cpu), r1, atol=1e-6)
        end
    end
end
