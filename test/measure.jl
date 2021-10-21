using Test, RydbergEmulator
using Random
using Graphs
using LinearAlgebra
using BitBasis

if !isdefined(@__MODULE__, :test_graph)
    include("utils.jl")
end

@testset "test measurement" begin
    Random.seed!(8)
    params = randn(10)
    ts = params[1:2:end]
    ϕs = params[2:2:end]
    hs = simple_rydberg.(nv(test_graph), ϕs)

    # prepare a zero state
    r = RydbergEmulator.zero_state(5, test_subspace)
    sample1 = measure!(r)
    @test sample1 == zero(BitStr64{5})
    emulate!(r, ts, hs)
    Random.seed!(5)
    # 1. sampling
    samples = measure(r; nshots=10000)
    @test samples isa Vector{<:BitStr}
    expected_sampling = samples .|> count_vertices |> mean
    # 2. exact
    expected_exact = r |> mean_rydberg
    @test isapprox(expected_exact, expected_sampling; rtol=1e-1)

    ####### gibbs
    Random.seed!(5)
    # 1. sampling
    expected_sampling = r |> measure(nshots=10000) |> gibbs_loss(0.5)
    # 2. exact
    expected_exact = r |> gibbs_loss(0.5)
    @test isapprox(expected_exact, expected_sampling; rtol=1e-1)
end
