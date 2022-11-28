using Test

if "docstring" in ARGS
    # include("docstrings.jl")
    exit()
end

@testset "problem" begin
    include("problem.jl")
end

@testset "forward diff" begin
    include("forwarddiff.jl")
end
