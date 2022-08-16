using Test
using BloqadeExpr
using BloqadeExpr: default_unit
using Unitful: kHz, µHz, Hz, MHz, µm, m, mm, cm

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
    @test h.rydberg_term.C ≈ 0.1092
end

unitless_values = [5, 0:0.1:1, [1,2,3]]
unit_values = [5m, 0cm:0.1cm:1cm, [1mm,2mm,3mm]]
conversion_exponent = [6, 4, 3]

@testset "convert $(typeof(unitful))" for (unitless, unitful, exponent) in zip(unitless_values, unit_values, conversion_exponent)
    @test default_unit(μm, unitful) == unitless*10^exponent
end
