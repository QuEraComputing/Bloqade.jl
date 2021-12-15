using Test
using RydbergEmulator
using LinearAlgebra
using Random
using BitBasis
using Yao
using Graphs
using RydbergEmulator: add_vertices!, add_random_vertices

if !isdefined(@__MODULE__, :test_graph)
    include("utils.jl")
end

@testset "loss functions" begin
    constraint_r = RydbergEmulator.zero_state(test_subspace)
    fullspace_r = Yao.zero_state(5)

    for loss_fn in [mean_rydberg, x->gibbs_loss(x, 0.3)]
        @test loss_fn(constraint_r) == 0.0
        @test loss_fn(fullspace_r) == 0.0
        @test loss_fn(measure(fullspace_r; nshots=10)) == 0.0
        @test loss_fn(measure(constraint_r; nshots=10)) == 0.0
    end
end

# generate random atom positions
atoms = RydAtom.([(0.0, 1.0), (1.0, 0.), (2.0, 0.0),
(1.0, 1.0), (1.0, 2.0), (2.0, 2.0)])
graph = unit_disk_graph(atoms, 1.5)
config = [1, 1, 1, 0, 1, 1]

@test RydbergEmulator.num_mis_violation(config, graph, 1) == 2
@test RydbergEmulator.num_mis_violation(config, graph, 2) == 2
@test RydbergEmulator.num_mis_violation(config, graph, 3) == 1
@test RydbergEmulator.num_mis_violation(config, graph, 4) == 0
@test RydbergEmulator.num_mis_violation(config, graph, 5) == 2
@test RydbergEmulator.num_mis_violation(config, graph, 6) == 1

@test !is_independent_set(config, graph)
to_independent_set!(config, graph)
@test is_independent_set(config, graph)


# config = bitarray(36, length(atoms))
# RydbergEmulator.to_independent_set!(config, graph)

# space = independent_set_subspace(graph)
# raw_state = zeros(ComplexF64, length(space))
# raw_state[space.map[packbits(config)]] = 1.0
# r = RydbergReg(length(atoms), raw_state, space)

# TODO: add an violation test

@testset "mis probabilities" begin
    raw_state = normalize!(rand(ComplexF64, length(space)))
    r = RydbergReg(raw_state, space)
    @test sum(independent_set_probabilities(r, graph)) ≈ 1
    @test sum(independent_set_probabilities(mis_postprocessing(graph), r, graph)) ≈ 1
end

@testset "RealLayout MIS functions" begin
    space = Subspace(10, sort!(randperm(1<<10)[1:30]))
    cr = rand_state(space)
    rr = RydbergReg{RealLayout}(cr)
    @test mean_rydberg(cr) ≈ mean_rydberg(rr)
end

@testset "mis_postprocessing" begin
    config = [0, 0, 0, 0, 0]
    @test add_vertices!(config, test_graph, 1:5) == [1,0,1,0,1]
    @test add_random_vertices(config, test_graph) == [1,0,1,0,1]
    @test count_vertices(mis_postprocessing(0, test_graph)) > 0        
end

@testset "SubspaceMap" begin
    atoms = square_lattice(10, 0.8)
    graph = unit_disk_graph(atoms, 1.5)
    space = independent_set_subspace(graph)
    reg = rand_state(space)
    Random.seed!(1234)
    l1 = mean_rydberg(mis_postprocessing(graph), reg)
    Random.seed!(1234)
    l2 = mean_rydberg(SubspaceMap(mis_postprocessing(graph), space), reg)
    @test l1 ≈ l2
end
