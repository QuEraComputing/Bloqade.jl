using Yao
using EaRydODE
using OrdinaryDiffEq
using Test

atoms = square_lattice(5, 0.8)
space = blockade_subspace(atoms, 1.5)

@testset "h=$name" for (name, h) in [
    "x+z" => XTerm(5, 1.0) + ZTerm(5, sin),
    "rydberg" => rydberg_h(atoms;Δ=sin, Ω=cos, C=2π * 109),
]

    dt = 1e-5
    @testset "subspace" begin
        ref = zero_state(space)
        discrete = KrylovEvolution(ref, map(_->dt, 0.0:dt:0.2), map(h, 0.0:dt:0.2))
        emulate!(discrete)

        reg = zero_state(space, ComplexLayout())
        continuous = ODEEvolution(reg, 0.2, h)
        emulate!(continuous)
        @test reg ≈ ref atol=1e-4

        reg = zero_state(space, RealLayout())
        emulate!(ODEEvolution(reg, 0.2, h))
        @test reg ≈ ref atol=1e-4
    end

    @testset "fullspace" begin
        ref = zero_state(5)
        emulate!(KrylovEvolution(ref, map(_->dt, 0.0:dt:0.2), map(h, 0.0:dt:0.2)))

        reg = zero_state(5)
        emulate!(ODEEvolution(reg, 0.2, h))
        @test ref ≈ reg atol=1e-4
    end
end

@testset "assertion" begin
    @test_throws ErrorException ODEEvolution(zero_state(9), 0.1, XTerm(5, 1.0))
end
