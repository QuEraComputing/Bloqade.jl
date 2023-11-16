using Test
using BloqadeExpr
using YaoArrayRegister
using BloqadeWaveforms
using BloqadeKrylov
using BloqadeLattices
using BloqadeExpr: Hamiltonian
using Yao


@testset "cfet42_sinX" begin
    
    atoms = generate_sites(ChainLattice(), 1, scale = 1)
    wf = Waveform(t->2.2*2π*sin(2π*t), duration = 1.3);
    h = rydberg_h(atoms; Ω = wf)

    ## do magnus4
    reg = zero_state(length(atoms))
    clocks = collect(0:1e-3:1.3)
    prob = CFETEvolution(reg, clocks, h, CFET4_2())
    show(stdout, MIME"text/plain"(), prob)
    emulate!(prob)

    @test prob.reg.state ≈ solution(1.3)

    reg = zero_state(length(atoms))
    prob = CFETEvolution(reg, clocks, h, CFET4_2())
    for info in prob
        @test info.reg.state ≈ solution(info.clock)
    end
    
end


@testset "cfet42_sinX expm multiply backend" begin
    
    atoms = generate_sites(ChainLattice(), 1, scale = 1)
    wf = Waveform(t->2.2*2π*sin(2π*t), duration = 1.3);
    h = rydberg_h(atoms; Ω = wf)

    ## do magnus4
    reg = zero_state(length(atoms))
    clocks = collect(0:1e-3:1.3)
    prob = CFETEvolution(reg, clocks, h, CFET4_2(), expmv_backend = expm_multiply!)
    show(stdout, MIME"text/plain"(), prob)
    #@test_throws ArgumentError KrylovEvolution(reg, [-0.1, 0.1], h)
    emulate!(prob)

    @test prob.reg.state ≈ solution(1.3)

    reg = zero_state(length(atoms))
    prob = CFETEvolution(reg, clocks, h, CFET4_2(), expmv_backend = expm_multiply!)
    for info in prob
        @test info.reg.state ≈ solution(info.clock)
    end

end
