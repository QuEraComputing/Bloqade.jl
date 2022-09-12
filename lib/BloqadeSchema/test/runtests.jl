using BloqadeSchema
using Configurations
using Test



@testset "parse" begin
    include("parse.jl")
end

@testset "transform" begin
    include("transform.jl")
end

# @testset "execute" begin
#     include("execute.jl")
# end

# @testset "serialize" begin
#     include("serialize.jl")
# end
