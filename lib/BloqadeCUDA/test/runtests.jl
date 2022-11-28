using Test

if "docstring" in ARGS
    # include("docstrings.jl")
    exit()
end

@testset "krylov" begin
    include("krylov.jl")
end

@testset "ode" begin
    include("ode.jl")
end
