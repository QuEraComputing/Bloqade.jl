using Test
using BloqadeExpr
using YaoBlocks
using BloqadeExpr: emit_dynamic_terms, emit_lowered, Hamiltonian, to_matrix, XPhase, PuPhase, PdPhase

@testset "emit_lowered" begin
    @test emit_lowered(SumOfX(5, 0.1)) == 0.1 * sum(put(5, i => X) for i in 1:5)
    @test emit_lowered(SumOfN(5, 0.1)) == 0.1 * sum(put(5, i => ConstGate.P1) for i in 1:5)
    @test emit_lowered(SumOfZ(5, 0.1)) == 0.1 * sum(put(5, i => Z) for i in 1:5)
    @test emit_lowered(SumOfXPhase(5, 0.1, 0.2)) == sum(0.1 * put(5, i => XPhase(0.2)) for i in 1:5)
end

atoms = [(1, 1), (1, 2), (1, 3)]
params = [nothing, 1.0, 2.0, [0.1 for _ in 1:3], sin, [sin for _ in 1:3], cos]
@testset "emit_dynamic_terms" begin
    @testset "Ω=$Ω" for Ω in params
        @testset "Δ=$Δ" for Δ in params
            @testset "ϕ=$ϕ" for ϕ in params
                h = rydberg_h(atoms; Ω, Δ, ϕ)
                H = Hamiltonian(Float64, h)
                @test to_matrix(H(0.1)) ≈ mat(Float64, h |> attime(0.1))
            end
        end
    end
end
