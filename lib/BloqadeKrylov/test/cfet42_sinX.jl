using Test
using BloqadeExpr
using YaoArrayRegister
using BloqadeWaveforms
using BloqadeKrylov
using BloqadeLattices
using BloqadeExpr: Hamiltonian
using BloqadeODE
using Yao

@testset "cfet42_sinX" begin

    atoms = generate_sites(ChainLattice(), 1, scale = 1)
    wf = Waveform(t->2.2*2π*sin(2π*t), duration = 1.3);
    h = rydberg_h(atoms; Ω = wf)

    ## do magnus4
    reg = zero_state(length(atoms))
    clocks = collect(0:1e-3:1.3)
    prob = CFET42Evolution(reg, clocks, h)
    show(stdout, MIME"text/plain"(), prob)
    #@test_throws ArgumentError KrylovEvolution(reg, [-0.1, 0.1], h)
    emulate!(prob)

    
    #benchmark against ODE solver:
    odereg = zero_state(length(atoms))
    ODEprob = SchrodingerProblem(odereg,1.3,h)
    show(stdout, MIME"text/plain"(), ODEprob)
    emulate!(ODEprob)

    @test prob.reg.state ≈ ODEprob.reg.state 

    prob = CFET42Evolution(reg, clocks, h)
    for info in prob
        @test info.clock == clocks[info.step]
    end
end
