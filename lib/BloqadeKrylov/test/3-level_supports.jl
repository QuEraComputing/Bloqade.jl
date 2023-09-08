using Test
using BloqadeExpr
using BloqadeKrylov
using Yao

@testset "The Levine-Pichler gate" begin
    nsites = 2
    V = 1e3
    Ω_r = 1.0
    τ = 4.29268/Ω_r
    Δ_r = 0.377371*Ω_r
    ϕ_r = 3.90242  # sign flipped
    atoms = [(0.0, 0.0), (0.0, (862690*2pi/V)^(1/6))]

    reg = zero_state(nsites; nlevel = 3)
    h1 = rydberg_h_3(atoms; Ω_hf = 1)
    prob1 = KrylovEvolution(reg, 0:1e-4:1/2*pi, h1)
    emulate!(prob1)
    goal1 = [0.5; -0.5im; 0; -0.5im; -0.5; 0; 0; 0; 0]  # ((|0> - i|1>)/√2)^{⊗2}
    @test isapprox(state(reg), goal1; atol = 1e-3)

    h2 = rydberg_h_3(atoms; Ω_r = Ω_r, Δ_r = Δ_r)
    prob2 = KrylovEvolution(reg, 0:1e-4:τ, h2)
    emulate!(prob2)
    goal2 = exp(-im*τ*Matrix(h2)) * goal1
    @test isapprox(state(reg), goal2; atol = 1e-3)

    h3 = rydberg_h_3(atoms; Ω_r = Ω_r, Δ_r = Δ_r, ϕ_r = ϕ_r)
    BloqadeExpr.mat(h3)
    prob3 = KrylovEvolution(reg, 0:1e-4:τ, h3)
    emulate!(prob3)
    goal3 = exp(-im*τ*Matrix(h3)) * goal2
    @test isapprox(state(reg), goal3; atol = 1e-3)
end

@testset "|11> -> -i/√2*(|1r> + |r1>)" begin
    nsites = 2
    V = 1e5
    atoms = [(0.0, 0.0), (0.0, (862690*2pi/V)^(1/6))]
    reg = product_state(dit"11;3")
    h = rydberg_h_3(atoms; Ω_r = 1)

    goal = exp(-im*pi/sqrt(2)*Matrix(h)) * state(reg)
    prob = KrylovEvolution(reg, 0:1e-2:pi/sqrt(2), h)
    emulate!(prob)
    @test isapprox(state(reg), goal; atol = 1e-2)
end
