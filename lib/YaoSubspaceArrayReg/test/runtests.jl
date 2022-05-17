using Test

@testset "register" begin
    include("type.jl")
end

@testset "instruct" begin
    include("instruct.jl")
end

@testset "measure" begin
    include("measure.jl")
end
