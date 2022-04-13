using Test
using BloqadeExpr
using YaoBlocks
using BloqadeExpr: emit_dynamic_terms, emit_lowered, Hamiltonian

@testset "emit_lowered" begin
    @test emit_lowered(SumOfX(5, 0.1)) == 0.1 * sum(put(5, i=>X) for i in 1:5)
    @test emit_lowered(SumOfN(5, 0.1)) == 0.1 * sum(put(5, i=>ConstGate.P1) for i in 1:5)
    @test emit_lowered(SumOfZ(5, 0.1)) == 0.1 * sum(put(5, i=>Z) for i in 1:5)
    @test emit_lowered(SumOfXPhase(5, 0.1, 0.2)) == sum(0.1 * put(5, i=>XPhase(0.2)) for i in 1:5)
end

@testset "emit_dynamic_terms" begin
    atoms = [(1, 1), (1, 2), (1, 3), (1, 4), (1, 5)]
    h = rydberg_h(atoms; Ω=1.0, Δ=sin)
    H = Hamiltonian(Float64, h)
    @test sum(zip(H.fs, H.ts)) do (f, h)
        f(0.1) * h
    end ≈ mat(Float64, h|>attime(0.1))        
end
