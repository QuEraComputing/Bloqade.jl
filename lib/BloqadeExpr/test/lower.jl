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
@testset "emit_dynamic_terms (2-level)" begin
    @testset "Ω=$Ω" for Ω in params
        @testset "Δ=$Δ" for Δ in params
            @testset "ϕ=$ϕ" for ϕ in params
                h = rydberg_h(atoms; Ω, Δ, ϕ)
                H = Hamiltonian(ComplexF64, h)
                @test to_matrix(H(0.1)) ≈ mat(ComplexF64, h |> attime(0.1))
            end
        end
    end
end;

nsites = 3
atoms = [(1, i) for i = 1:nsites]
params = [nothing, 1.0, [0.1 for _ in 1:nsites], sin, [sin for _ in 1:nsites]]
@testset "emit_dynamic_terms (3-level)" begin
    @testset "Ω_hf=$Ω_hf, Δ_hf=$Δ_hf, ϕ_hf=$ϕ_hf" for Ω_hf in params[2:2], Δ_hf in params[3:4], ϕ_hf in params[5:5]
        @testset "Ω_r=$Ω_r, Δ_r=$Δ_r, ϕ_r=$ϕ_r" for Ω_r in params[2:end], Δ_r in params, ϕ_r in params
            h = rydberg_h_3(atoms; Ω_hf, Δ_hf, ϕ_hf, Ω_r, Δ_r, ϕ_r)
            H = Hamiltonian(ComplexF64, h)
            @test to_matrix(H(0.1)) ≈ mat(ComplexF64, h |> attime(0.1))
        end
    end
    @testset "SumOfZ_01, SumOfZ_1r" for a in params, b in params
        h = SumOfZ_01(nsites, isnothing(a) ? 0.0 : a) + SumOfZ_1r(nsites, isnothing(b) ? 0.0 : b)
        H = Hamiltonian(ComplexF64, h)
        @test to_matrix(H(0.1)) ≈ mat(ComplexF64, h |> attime(0.1))
    end
end;
