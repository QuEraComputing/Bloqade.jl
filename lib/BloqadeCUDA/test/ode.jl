using CUDA
using Adapt
using BloqadeODE
using BloqadeCUDA
using BloqadeWaveforms
using BloqadeLattices
using YaoArrayRegister
CUDA.allowscalar(false)

@testset "ODE solver" begin
    nsites = 5
    total_time = 3.0
    Ω_max = 2π * 4
    Ω = piecewise_linear(clocks = [0.0, 0.1, 2.1, 2.2, total_time], values = [0.0, Ω_max, Ω_max, 0, 0])
    U1 = -2π * 10
    U2 = 2π * 10
    Δ = piecewise_linear(clocks = [0.0, 0.6, 2.1, total_time], values = [U1, U1, U2, U2])
    atoms = generate_sites(ChainLattice(), nsites, scale = 5.72)
    h = rydberg_h(atoms; Δ, Ω)
    reg = zero_state(nsites)
    prob = SchrodingerProblem(reg, total_time, h)
    emulate!(prob)
    d_prob = adapt(CuArray, prob)
    emulate!(d_prob)

    @test Array(d_prob.reg.state) ≈ prob.reg.state
end
