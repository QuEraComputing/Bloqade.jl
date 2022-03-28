using Test
using EaRydExpr

@testset "rydberg_h" begin
    positions = [(1, 2), (2, 3)]
    @test rydberg_h(positions; Ω=1.0) ==
        RydInteract(;atoms=positions) + SumOfX(;nsites=2, Ω=0.5)

    @test rydberg_h(positions; Ω=1.0, Δ=0.2) ==
        RydInteract(;atoms=positions) + SumOfX(;nsites=2, Ω=0.5) -
            SumOfN(;nsites=2, Δ=0.2)
end
