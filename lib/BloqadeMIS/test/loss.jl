using Test
using BloqadeMIS
using LinearAlgebra
using Random
using BitBasis
using YaoSubspaceArrayReg
using Graphs
using BloqadeMIS: add_vertices!, add_random_vertices
using FromFile
using BloqadeLattices
using Statistics

Random.seed!(1234)

@from "utils.jl" import test_subspace, test_graph

@testset "loss functions" begin
    constraint_r = zero_state(test_subspace)
    fullspace_r = zero_state(5)

    for loss_fn in [rydberg_density_sum, x->gibbs_loss(x, 0.3)]
        @test loss_fn(constraint_r) == 0.0
        @test loss_fn(fullspace_r) == 0.0
        @test loss_fn(measure(fullspace_r; nshots=10)) == 0.0
        @test loss_fn(measure(constraint_r; nshots=10)) == 0.0
    end
end

# generate random atom positions
atoms = [(0.0, 1.0), (1.0, 0.), (2.0, 0.0), (1.0, 1.0), (1.0, 2.0), (2.0, 2.0)]
graph = unit_disk_graph(atoms, 1.5)
config = [1, 1, 1, 0, 1, 1]

@test BloqadeMIS.num_mis_violation(config, graph, 1) == 2
@test BloqadeMIS.num_mis_violation(config, graph, 2) == 2
@test BloqadeMIS.num_mis_violation(config, graph, 3) == 1
@test BloqadeMIS.num_mis_violation(config, graph, 4) == 0
@test BloqadeMIS.num_mis_violation(config, graph, 5) == 2
@test BloqadeMIS.num_mis_violation(config, graph, 6) == 1

@test !is_independent_set(config, graph)
to_independent_set!(config, graph)
@test is_independent_set(config, graph)

@testset "mis probabilities" begin
    raw_state = normalize!(rand(ComplexF64, length(test_subspace)))
    r = SubspaceArrayReg(raw_state, test_subspace)
    @test sum(independent_set_probabilities(r, graph)) ≈ 1
    @test sum(independent_set_probabilities(mis_postprocessing(graph), r, graph)) ≈ 1
end

@testset "mis_postprocessing" begin
    config = [0, 0, 0, 0, 0]
    @test add_vertices!(config, test_graph, 1:5) == [1,0,1,0,1]
    @test add_random_vertices(config, test_graph) == [1,0,1,0,1]
    @test count_vertices(mis_postprocessing(0, test_graph)) > 0        
end


@testset "SubspaceMap" begin
    atoms = generate_sites(SquareLattice(), 4, 4) |> random_dropout(0.2)
    graph = unit_disk_graph(atoms, 1.5)
    space = independent_set_subspace(graph)
    reg = rand_state(space)
    Random.seed!(1234)
    l1 = rydberg_density_sum(mis_postprocessing(graph), reg)
    Random.seed!(1234)
    l2 = rydberg_density_sum(SubspaceMap(mis_postprocessing(graph), space), reg)
    @test l1 ≈ l2
end

@testset "exact/sample based loss function" begin
    r = rand_state(5)
    samples = measure(r; nshots=10000)
    expected_sampling = samples .|> count_vertices |> mean
    # 2. exact
    expected_exact = r |> rydberg_density_sum
    @test isapprox(expected_exact, expected_sampling; rtol=1e-1)

    ####### gibbs
    Random.seed!(5)
    # 1. sampling
    expected_sampling = r |> measure(nshots=10000) |> gibbs_loss(0.5)
    # 2. exact
    expected_exact = r |> gibbs_loss(0.5)
    @test isapprox(expected_exact, expected_sampling; rtol=1e-1)
end
