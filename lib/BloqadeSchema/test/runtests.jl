using BloqadeSchema
using Configurations
using Test

@testset "execute" begin
    include("execute.jl")
end

@testset "parse" begin
    include("parse.jl")
end

# @testset "serialize" begin
#     include("serialize.jl")
# end
