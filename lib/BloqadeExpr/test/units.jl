using Test
using BloqadeExpr
using Unitful: kHz, µm, m, mm, cm
using BloqadeExpr:default_unit

@testset "term units" begin
    @test SumOfXPhase(5, 2kHz, 1.0).Ω ≈ 0.002
    @test SumOfX(5, 2kHz).Ω ≈ 0.002
    @test SumOfZ(5, 2kHz).Δ ≈ 0.002
    @test SumOfN(5, 2kHz).Δ ≈ 0.002

    h = RydInteract(atoms = [(1,), (2,), (3,), (4)], C = 109.2kHz * µm^6)
    @test h.C ≈ 0.1092

    h = rydberg_h([(1,), (2,)], C = 109.2kHz * µm^6)
    @test h.C ≈ 0.1092
end

unitless_values = [5, 0:0.1:1, [1,2,3]]
unit_values = [5m, 0cm:0.1cm:1cm, [1mm,2mm,3mm]]
conversion_exponent = [6, 4, 3]

#test for a value::Quantity, range::AbstractRange and vector::AbstractArray{S} where S<:Quantity
@testset "basic unit conversion tests" for i in 1:3
    @test default_unit(μm, unit_values[i]) == unitless_values[i]*10^conversion_exponent[i]
end
