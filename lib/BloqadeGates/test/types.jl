using BloqadeGates
using BloqadeExpr, BloqadeODE
using Yao
using Test
using BloqadeGates: two_level_indices

atoms = [(0.0, 0.0), (10.0, 0.0)]
Ω, ϕ, Δ = rand(3)
rh = rydberg_h(atoms; Ω, ϕ, Δ)
rh3 = rydberg_h_3(atoms; Ω_hf = Ω, ϕ_hf = ϕ, Δ_hf = Δ, Ω_r = Ω, ϕ_r = ϕ, Δ_r = Δ)
@testset "mat" begin
    # different starting time
    @test mat(RydbergPulse(rh, 1.0)) ≈ mat(RydbergPulse(rh, 1.0, 2.0)) 
    @test mat(RydbergPulse(rh3, 1.0)) ≈ mat(RydbergPulse(rh3, 1.0, 2.0))

    # different complex type
    @test mat(RydbergPulse(rh, 1.0)) ≈ mat(ComplexF32, RydbergPulse(rh, 1.0))
    @test mat(RydbergPulse(rh3, 1.0)) ≈ mat(ComplexF32, RydbergPulse(rh3, 1.0))

    # time-dependent Hamiltonian
    fid = operator_fidelity(RydbergPulse(rydberg_h(atoms; Ω = one), pi; backend = SchrodingerProblem), 
        RydbergPulse(rydberg_h(atoms; Ω = 1.0), pi; backend = SchrodingerProblem))
    @test 1-1e-4 < fid < 1+1e-4
    ids = two_level_indices(2)
    p_t = RydbergPulse(rydberg_h_3(atoms; Ω_hf = one), pi; step = 1e-3)
    p_1 = RydbergPulse(rydberg_h_3(atoms; Ω_hf = 1.0), pi; step = 1e-3)
    @test operator_fidelity(matblock(mat(p_t)[ids, ids]), matblock(mat(p_1)[ids, ids])) > 1-1e-6
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