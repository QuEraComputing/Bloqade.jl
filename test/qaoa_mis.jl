using Test
using RydbergEmulator
using Yao
using LightGraphs

if !isdefined(@__MODULE__, :test_graph)
    include("utils.jl")
end

@testset "loss functions" begin
    constraint_r = RydbergEmulator.zero_state(nv(test_graph), test_subspace)
    fullspace_r = Yao.zero_state(5)

    for loss_fn in [mean_rydberg, x->gibbs_loss(x, 0.3)]
        @test loss_fn(constraint_r) == 0.0
        @test loss_fn(fullspace_r) == 0.0
        @test loss_fn(measure(fullspace_r; nshots=10)) == 0.0
        @test loss_fn(measure(constraint_r; nshots=10)) == 0.0
    end    
end
