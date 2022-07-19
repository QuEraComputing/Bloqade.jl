using Test
using BloqadeExpr
using BloqadeExpr: default_unit
using Unitful: kHz, µHz, Hz, MHz, µm

@testset "term units" begin
    @test SumOfXPhase(5, 2kHz, 1.0).Ω ≈ 0.002
    @test SumOfX(5, 2kHz).Ω ≈ 0.002
    @test SumOfZ(5, 2kHz).Δ ≈ 0.002
    @test SumOfN(5, 2kHz).Δ ≈ 0.002
    @test SumOfXPhase(5, 2kHz, 2.0).ϕ ≈ 2.0

    @test SumOfX(3, [2μHz,2Hz,2MHz]).Ω ≈ [2*1e-12,2*1e-6,2]

    h = RydInteract(atoms = [(1,), (2,), (3,), (4)], C = 109.2kHz * µm^6)
    @test h.C ≈ 0.1092
    h = rydberg_h([(1,), (2,)], C = 109.2kHz * µm^6)
    @test h.C ≈ 0.1092
end

