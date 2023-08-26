using BloqadeNoisy
using BloqadeExpr
using BloqadeWaveforms
using Kronecker
using StatsBase
using SciMLBase
using YaoArrayRegister
using YaoBlocks
using LinearAlgebra

@testset "problem interface" begin
    reg = zero_state(1)
    h = rydberg_h([(0, 0)], Ω=2π, Δ=0)
    save_times = LinRange(0.0, 1.0, 10)
    ns = @test_nowarn NoisySchrodingerProblem(reg, save_times, h, Aquila())
end

@testset "emulation" begin
    save_times = LinRange(0.0, 1.0, 10)
    reg = zero_state(1)
    ns = NoisySchrodingerProblem(
        reg, save_times, rydberg_h([(0, 0)], Ω=2π, Δ=0), Aquila()
    )
    sim = @test_nowarn emulate_noisy(ns, 1)
    @test length(sim) == length(save_times)
    @test length(first(sim)) == length(reg.state)
end

@testset "Hamiltonian sampling" begin
    h = rydberg_h(
        [(0.0,)];
        Ω = Waveform(sin, 4),
        Δ = piecewise_constant(
            clocks = [0, 4],
            values = [2π],
        ),
        ϕ = piecewise_constant(
            clocks = [0, 4],
            values = [0]
        )
    )
    sampler = Aquila().coherent_noise(h)
    @test_nowarn solve(SchrodingerProblem(zero_state(1), 4.0, sampler()), DP8())
end

@testset "noisy measurements" begin
    state = [[.5]; [.5]]
    @test length(measure_noisy(Aquila(), state; nshots=10)) == 10
end

@testset "expectation values" begin
    ns = NoisySchrodingerProblem(
        zero_state(1), [0.0, 4.0], rydberg_h([(0, 0)], Ω=2π, Δ=0), Aquila()
    )
    sim = @test_nowarn emulate_noisy(ns, 1, mat.([X, Z]))
    @test sim[1][1] == 0.0
    @test sim[2][1] == 1.0
end

@testset "noisy expectation values" begin
    ns = NoisySchrodingerProblem(
        zero_state(1), [0.0, 4.0], rydberg_h([(0, 0)], Ω=2π, Δ=0), Aquila()
    )
    p01 = Aquila().confusion_mat(1)[2, 1]
    sim = emulate_noisy(ns, 1, [mat(Z)]; readout_error=true)
    @test sim[1][1] == 1 - 2 * p01
    p = [[.5]; [.5]]
    p_ro = Aquila().confusion_mat(1) * p
    Zexpec = expectation_value_noisy(Aquila(), p, mat(Z))
    @test Zexpec == p_ro[1] - p_ro[2]
    Zexpec = expectation_value_noisy(Aquila(), p, mat(Z), errs = [.05, .05])
    e_t = Aquila().confusion_mat(1) * [[.05]; [.05]]
    @test Zexpec.propagated_err == norm(e_t)
    Z_shots = expectation_value_noisy(Aquila(), p, mat(Z), 1000)
    @test isapprox(Z_shots, Zexpec.expectation, atol = 4/sqrt(1000))
    @test expectation_value_noisy(Aquila(), p, mat(Z), 1000; errs = true).sample_err < .1
end

@testset "emulation arguments" begin
    h = rydberg_h([(0, 0), (8, 0), (18, 0)], Ω=15, Δ=0)
    reg = zero_state(3)
    save_times = LinRange(0.0, 1.0, 5)
    ns = NoisySchrodingerProblem(reg, save_times, h, Aquila())
    sim1 = @test_nowarn emulate_noisy(ns, 1, [mat(put(3, i => Z)) for i in 1:3];
        report_error=true
    )
    sim2 = @test_nowarn emulate_noisy(ns, 2, [mat(put(3, 1 => Z))];
        report_error=true, ensemble_algo=EnsembleThreads()
    )
    sim_ro = @test_nowarn emulate_noisy(ns, 2, [mat(put(3, 1 => Z))];
        readout_error=true, report_error=true
    )
    sim = @test_nowarn emulate_noisy(ns, 1, [mat(put(3, i => Z)) for i in 1:3];
        readout_error=true, report_error=true, shots=10
    )
end

@testset "ensemble analysis" begin
    ns = NoisySchrodingerProblem(
        zero_state(1), [0.0, 4.0], rydberg_h([(0, 0)], Ω=2π, Δ=0), Aquila()
    )
    sim = emulate_noisy(ns, 2, sol -> [abs(u[1])^2 for u in sol])
    @test values = simulation_series_mean(sim; index = 1)[end] ==
                   mean([s[end] for s in sim])
    @test simulation_series_err(sim; index=1, factor=3)[end] ==
          std([s[end][1] for s in sim]) * 3 / sqrt(2)
end

@testset "EnsembleProblem integration" begin
    ns = NoisySchrodingerProblem(
        zero_state(1), [0.0, 4.0], rydberg_h([(0, 0)], Ω=2π, Δ=0), Aquila()
    )
    @test_nowarn randomize(ns)
    ep = @test_nowarn EnsembleProblem(ns, 
        prob_func=(prob, i, repeat) -> randomize(ns)
    )
    @test_nowarn solve(ep, DP8(); trajectories=1)
end

@testset "trivial noise model" begin
    reg = zero_state(1)
    save_times = [0, 4]
    h = rydberg_h([(0, 0)]; Ω=2π)
    trivial_error_model = ErrorModel(
        n -> I,
        n -> [],
        h -> (() -> h)
    )
    ns = NoisySchrodingerProblem(reg, save_times, h, trivial_error_model)
    sim = emulate_noisy(ns, 10, sol -> sol[end])
    ψ = first(solve(SchrodingerProblem(reg, 4, h), DP8()).u)
    @test all(map(sim) do u; u ≈ ψ; end)
end

@testset "custom collapse operators" begin
    save_times = [0, 4]
    h = rydberg_h([(0, 0)]; Ω=2π)
    rate = 1 / 10
    c_ops = [sqrt(rate) * mat((X + im * Y) / 2)]
    @test_nowarn ns = NoisySchrodingerProblem(zero_state(1), save_times, h, c_ops)
end

@testset "custom error model" begin
    confusion_matrix(n) = kronecker([[[0.9 0.1]; [0.1 0.9]] for i in 1:n]...)
    bitflip_model(n) = [
        SparseMatrixCSC(sqrt(1 / 10) * mat(put(n, i => X))) for i in 1:n
    ]
    coherent_noise(h) = () -> (
        (atoms, ϕ, Ω, Δ) = get_rydberg_params(h);
        rydberg_h(atoms; Ω=Ω * (1 + 0.08 * randn()), Δ=Δ, ϕ=ϕ)
    )
    better_error_model = ErrorModel(
        confusion_matrix,
        bitflip_model,
        coherent_noise
    )
    ns = NoisySchrodingerProblem(zero_state(2), 0:1.0f-2:1, rydberg_h([(0,), (8)]; Ω=15), better_error_model)
    sim = @test_nowarn emulate_noisy(ns, 1, [mat(put(2, 1 => X)), mat(kron(X, X))])
end