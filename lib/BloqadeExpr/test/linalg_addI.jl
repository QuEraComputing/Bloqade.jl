using Test
using BloqadeExpr
using LinearAlgebra
using Random
using LuxurySparse
using BloqadeExpr: Hamiltonian

atoms = [(1, 1), (1, 2), (1, 3), (1, 4), (1, 5)]
h = rydberg_h(atoms; Ω = 1.0, Δ = sin)

hlist = Hamiltonian(Float64, h)


@testset "add_I for SumOfLinop" begin

    Ht = hlist(0.1)

    @test ishermitian(Ht)
    # these values are not mutated and can be re-used throughout the unit tests
    H1 = add_I(Ht, -0.5)
    H2 = add_I(Ht, -0.5+1im)
    H3 = add_I(Ht, 0.5im)

    @test ishermitian(H1) == true
    @test isskewhermitian(H1) == false

    @test ishermitian(H2) == false
    @test isskewhermitian(H2) == false

    @test ishermitian(H3) == false
    @test isskewhermitian(H3) == false


    Ha = 1.0im*Ht

    @test ishermitian(Ha) == false
    @test isskewhermitian(Ha) == true

    H1a = add_I(Ha, -0.5)
    H2a = add_I(Ha, -0.5+1im)
    H3a = add_I(Ha, 0.5im)   

    @test ishermitian(H1a) == false
    @test isskewhermitian(H1a) == false

    @test ishermitian(H2a) == false
    @test isskewhermitian(H2a) == false

    @test ishermitian(H3a) == false
    @test isskewhermitian(H3a) == true




    M = to_matrix(Ht)
    S = add_I(M, 1+3.)
    @test S == M + (1+3.)*LinearAlgebra.I(size(M,1))   



    
end