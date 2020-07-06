using Test, RydbergEmulator
using Random
using LightGraphs, LinearAlgebra
using BitBasis
using Yao

@testset "test measurement" begin
    Random.seed!(8)
    params = randn(10)
    ts = params[1:2:end]
    ϕs = params[2:2:end]
    hs = simple_rydberg.(nv(test_graph), ϕs)

    # prepare a zero state
    r = RydbergEmulator.zero_state(5, test_subspace)
    sample1 = Yao.measure!(r)
    @test sample1 == zero(BitStr64{5})
    emulate!(r, ts, hs)
    Random.seed!(5)
    # 1. sampling
    samples = measure(r; nshots=10000)
    @test samples isa Vector{<:BitStr}
    expected_sampling = samples .|> count_vertices |> mean
    # 2. exact
    expected_exact = r |> mean_nv
    @test isapprox(expected_exact, expected_sampling; rtol=1e-1)

    ####### softmis
    Random.seed!(5)
    # 1. sampling
    expected_sampling = r |> measure(nshots=10000) |> soft_misloss(0.5)
    # 2. exact
    expected_exact = r |> soft_misloss(0.5)
    @test isapprox(expected_exact, expected_sampling; rtol=1e-1)
end
