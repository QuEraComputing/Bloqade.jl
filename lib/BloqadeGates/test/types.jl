using BloqadeGates
using BloqadeExpr, BloqadeODE
using Yao
using Test

atoms = [(0.0, 0.0), (4.0, 0.0)]
Ω, ϕ, Δ = rand(3)
rh = rydberg_h(atoms; Ω, ϕ, Δ)
rh3 = rydberg_h_3(atoms; Ω_hf = Ω, ϕ_hf = ϕ, Δ_hf = Δ, Ω_r = Ω, ϕ_r = ϕ, Δ_r = Δ)
@testset "t_start, t_end" begin
    @test mat(RydbergPulse(rh, 1.0)) ≈ mat(RydbergPulse(rh, 1.0, 2.0))
    @test mat(RydbergPulse(rh3, 1.0)) ≈ mat(RydbergPulse(rh3, 1.0, 2.0))
end

@testset "ODE emulator" begin
    pulse_ode = RydbergPulse(rh, 1.0; backend = SchrodingerProblem)
    reg1 = rand_state(2)
    reg2 = copy(reg1)
    reg1 |> RydbergPulse(rh, 0.0, 1.0)
    reg2 |> pulse_ode
    normalize!(reg1)
    normalize!(reg2)
    @test 1-1e-4 < fidelity(reg1, reg2) < 1+1e-4
end

@testset "YaoAPI" begin
    @test nqudits(RydbergPulse(rh, 1.0)) == nqudits(RydbergPulse(rh, 1.0, 2.0)) == length(atoms)
    @test nqudits(RydbergPulse(rh3, 1.0)) == nqudits(RydbergPulse(rh3, 1.0, 2.0)) == length(atoms)
    @test nlevel(RydbergPulse(rh, 1.0)) == nlevel(RydbergPulse(rh, 1.0, 2.0)) == 2
    @test nlevel(RydbergPulse(rh3, 1.0)) == nlevel(RydbergPulse(rh3, 1.0, 2.0)) == 3
    @test occupied_locs(RydbergPulse(rh, 1.0)) == occupied_locs(RydbergPulse(rh, 1.0, 2.0)) == (1, 2)
    @test occupied_locs(RydbergPulse(rh3, 1.0)) == occupied_locs(RydbergPulse(rh3, 1.0, 2.0)) == (1, 2)
end

@testset "show" begin
    println(RydbergPulse(rh3, 1.0))
end