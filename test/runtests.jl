using Test
using Documenter
using Bloqade
using BloqadeLattices

doctest(BloqadeLattices; manual=false)

@testset "observables" begin
    include("observables.jl")    
end
