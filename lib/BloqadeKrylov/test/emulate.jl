using Test
using BloqadeExpr
using YaoArrayRegister
using BloqadeWaveforms
using BloqadeKrylov
using BloqadeLattices
using BloqadeExpr: Hamiltonian
using Yao

@testset "emulate" begin
    atoms = generate_sites(ChainLattice(), 4, scale = 3.2)
    clocks = [0.0, 0.5, 0.8, 1.1, 1.5, 1.6]
    wf = piecewise_constant(clocks = clocks, values = [0.0, 2.1, 2.1, 1.5, 0.0])
    h = rydberg_h(atoms; Ω = wf)
    reg = zero_state(length(atoms))
    prob = KrylovEvolution(reg, clocks, h)
    show(stdout, MIME"text/plain"(), prob)
    @test_throws ArgumentError KrylovEvolution(reg, [-0.1, 0.1], h)
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
