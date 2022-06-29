using Test
using BloqadeMIS

@testset "unit_disk_graph" begin
    include("unit_disk_graph.jl")
end

@testset "utils" begin
    include("utils.jl")
end

@testset "loss" begin
    include("loss.jl")
end
