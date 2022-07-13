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

#test for a value::Quantity
value = 5
value_m = 5m
@which default_unit(μm, value_m)
@test default_unit(μm, value_m) == value*10^6

#test for a range::AbstractRange
range = 0:0.1:1
range_cm = 0cm:0.1cm:1cm
@test default_unit(µm, range_cm) == range*10^4

#test for a vector::AbstractArray{S} where S<:Quantity
vector = [1,2,3]
vector_mm = [1mm,2mm,3mm]
@test default_unit(μm, vector_mm) == vector*10^3 



