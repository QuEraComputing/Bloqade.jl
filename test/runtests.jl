using Test
using Documenter
using Bloqade
using BloqadeExpr
using BloqadeLattices

# doctest(BloqadeExpr; manual=false)
# doctest(BloqadeLattices; manual = false)

if "docstring" in ARGS
    # include("docstrings.jl")
    exit()
end

@testset "observables" begin
    include("observables.jl")
end

@testset "plots" begin
    include("plots.jl")
end
