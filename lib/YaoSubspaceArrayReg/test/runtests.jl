using Test

if "docstring" in ARGS
    # include("docstrings.jl")
    exit()
end
@testset "register" begin
    include("type.jl")
end

@testset "instruct" begin
    include("instruct.jl")
end

@testset "measure" begin
    include("measure.jl")
end
