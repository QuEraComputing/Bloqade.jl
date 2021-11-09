using Yao
using ContinuousEmulator
using OrdinaryDiffEq
using Test

atoms = square_lattice(5, 0.8)
space = blockade_subspace(atoms, 1.5)

@testset "h=$name" for (name, h) in [
    "x+z" => XTerm(5, 1.0) + ZTerm(5, sin),
    "rydberg" => rydberg_h(atoms, sin, nothing, cos),
]

    dt = 1e-5
    @testset "subspace" begin
        ref = zero_state(space)
        emulate!(ref, map(_->dt, 0.0:dt:0.2), map(h, 0.0:dt:0.2))

        reg = zero_state(space, ComplexLayout())
        emulate!(reg, 0.2, h)
        @test reg ≈ ref atol=1e-4

        reg = zero_state(space, RealLayout())
        emulate!(reg, 0.2, h)
        @test reg ≈ ref atol=1e-4
    end

    @testset "fullspace" begin
        ref = zero_state(5)
        emulate!(ref, map(_->dt, 0.0:dt:0.2), map(h, 0.0:dt:0.2))

        reg = zero_state(5)
        emulate!(reg, 0.2, h)
        @test ref ≈ reg atol=1e-4
    end
end
