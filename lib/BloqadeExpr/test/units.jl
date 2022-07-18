using Test
using BloqadeExpr
using Unitful: kHz, µHz, Hz, µm, m

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

@testset "further units" begin
    @test SumOfXPhase(5, 2kHz, 2.0).ϕ ≈ 2.0 # Should test line 4, i.e. NoUnits version of default_unit
    @test default_unit(μm, (1.0:0.2:2.0)m) ≈ (1.0:0.2:2.0)* 1e6 # Should test lines 11-16, i.e. Range version of default unit. Note that Range does not seem to be used commonly in Hamiltonian terms (e.g. SumOfX, except for in Adiabatic Evolution)
    @test SumOfX(3, [2μHz,2Hz,2MHz]).Ω ≈ [2*1e-12,2*1e-6,2] # Should test lines 18-24, i.e. the Array version of default_unit
end
