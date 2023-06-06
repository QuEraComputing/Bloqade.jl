using Test
using BloqadeExpr
using YaoArrayRegister
using BloqadeWaveforms
using BloqadeKrylov
using BloqadeLattices
using BloqadeExpr: Hamiltonian
using Yao

@testset "krylov_sinX" begin

    atoms = generate_sites(ChainLattice(), 1, scale = 1)
    wf = Waveform(t->2.2*2π*sin(2π*t), duration = 1.3);
    h = rydberg_h(atoms; Ω = wf)
    reg = zero_state(length(atoms))
    clocks = collect(0:1e-2:1.3)
    prob = KrylovEvolution(reg, clocks, h)
    show(stdout, MIME"text/plain"(), prob)
    #@test_throws ArgumentError KrylovEvolution(reg, [-0.1, 0.1], h)
    emulate!(prob)

    hs = map(clocks[1:end-1]) do t
        return Matrix(h |> attime(t))
    end

    r = zero_state(length(atoms))
    state = statevec(r)
    durations = diff(clocks)
    for (t, h) in zip(durations, hs)
        state = exp(-im * t * h) * state
    end

    @test state ≈ prob.reg.state

    prob = KrylovEvolution(reg, clocks, h)
    for info in prob
        @test info.clock == clocks[info.step]
    end
end
