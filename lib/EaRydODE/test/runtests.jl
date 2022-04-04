using Test

@testset "problem" begin
    include("problem.jl")    
end

@testset "forward diff" begin
    include("forwarddiff.jl")
end
