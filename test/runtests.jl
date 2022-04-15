using Test
using Documenter
using Bloqade
using BloqadeExpr
using BloqadeLattices

doctest(BloqadeExpr; manual=false)
doctest(BloqadeLattices; manual=false)

@testset "observables" begin
    include("observables.jl")    
end
