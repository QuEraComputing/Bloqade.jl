using BloqadeSchema
using Test
using Configurations

@testset "execute" begin
    include("execute.jl")
end

@testset "serialize" begin
    include("serialize.jl")
end
