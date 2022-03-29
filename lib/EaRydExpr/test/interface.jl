using Test
using EaRydExpr
using YaoBlocks.Optimise

@testset "rydberg_h" begin
    positions = [(1, 2), (2, 3)]
    h = RydInteract(;atoms=positions) + SumOfX(;nsites=2, Ω=0.5)
    @test rydberg_h(positions; Ω=1.0) == Optimise.simplify(h)

    h = RydInteract(;atoms=positions) + SumOfX(;nsites=2, Ω=0.5) -
        SumOfN(;nsites=2, Δ=0.2)
    @test rydberg_h(positions; Ω=1.0, Δ=0.2) == Optimise.simplify(h)

    h = RydInteract(;atoms=positions) + SumOfXPhase(;nsites=2, Ω=0.5, ϕ=0.1) -
        SumOfN(;nsites=2, Δ=0.2)
    @test rydberg_h(positions; Ω=1.0, ϕ=0.1, Δ=0.2) == Optimise.simplify(h)

    h = RydInteract(;atoms=positions) + SumOfX(;nsites=2, Ω=[2.0, 1.0])
    @test rydberg_h(positions; Ω=[4.0, 2.0]) == Optimise.simplify(h)    
end


@testset "attime" begin
    h1 = RydInteract(;atoms=positions) + SumOfX(;nsites=2, Ω=sin)
    h2 = RydInteract(;atoms=positions) + SumOfX(;nsites=2, Ω=sin(0.1))
    @test h1 |> attime(0.1) == h2

    h1 = RydInteract(;atoms=positions) + SumOfX(;nsites=2, Ω=sin) -
        SumOfN(;nsites=2, Δ=cos)
    h2 = RydInteract(;atoms=positions) + SumOfX(;nsites=2, Ω=sin(0.1)) -
        SumOfN(;nsites=2, Δ=cos(0.1))
    @test h1 |> attime(0.1) == h2

    h1 = RydInteract(;atoms=positions) + SumOfXPhase(;nsites=2, Ω=0.5, ϕ=cos) -
        SumOfN(;nsites=2, Δ=0.2)
    h2 = RydInteract(;atoms=positions) + SumOfXPhase(;nsites=2, Ω=0.5, ϕ=cos(0.1)) -
        SumOfN(;nsites=2, Δ=0.2)
    @test h1 |> attime(0.1) == h2

    h1 = RydInteract(;atoms=positions) + SumOfX(;nsites=2, Ω=[sin, cos])
    h2 = RydInteract(;atoms=positions) + SumOfX(;nsites=2, Ω=[sin(0.1), cos(0.1)])
    @test h1 |> attime(0.1) == h2
end
