using Test
using BloqadeExpr
using LinearAlgebra
using Random
using LuxurySparse
using BloqadeExpr: Hamiltonian




@testset "derivative for Hamiltonain" begin

    atoms = [(1, 1)]
    h = rydberg_h(atoms; Ω = 1.0, Δ = sin)
    hlist = Hamiltonian(Float64, h)

    h2 = rydberg_h(atoms; Ω = 0.0,  Δ = cos)
    htar = Hamiltonian(Float64, h2)

    println(htar.ts)
    println(hlist.ts)

    for t in 0:0.1:2
        src = derivative(hlist, t)
        tar = htar(t)

        @test to_matrix(src) == to_matrix(tar)
    
    end


    
end