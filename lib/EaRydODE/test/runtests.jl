using Yao
using EaRydODE
using EaRydCore
using Test

@testset "ODE problem interface" begin
    reg = zero_state(5)
    tspan = (0, 1e-4)
    atoms = square_lattice(5, 0.8)
    h = rydberg_h(atoms; C=2π * 109.16, Ω=sin, ϕ=cos)
    prob = SchrodingerProblem(reg, tspan, h; dt=1e-5, progress=true, save_start=false)

    integrator = init(prob, Vern8())
    for (u, t) in tuples(integrator)
        @test u === prob.state
        @test vec(prob.reg.state) ≈ u
    end

    # integrator should be initialize properly
    integrator = init(prob, Vern8())
    @test integrator.u ≈ statevec(zero_state(5))
end

@testset "select specific time" begin
    reg = zero_state(5)
    tspan = (0, 0.2)
    atoms = square_lattice(5, 0.8)
    h = rydberg_h(atoms; C=2π * 109.16, Ω=sin, ϕ=cos)
    prob = SchrodingerProblem(reg, tspan, h; dt=1e-5, progress=true, save_start=false)
    integrator = init(prob, Vern8())
    for (u,t) in TimeChoiceIterator(integrator, [0.05, 0.1])
        @test t in [0.05, 0.1]
    end
end

@testset "ODE evolution results" begin
    atoms = square_lattice(5, 0.8)
    space = blockade_subspace(atoms, 1.5)

    @testset "h=$name" for (name, h) in [
        "x+z" => XTerm(5, 1.0) + NTerm(5, sin),
        "rydberg" => rydberg_h(atoms;Δ=sin, Ω=cos, C=2π * 109),
    ]

        dt = 1e-5
        @testset "subspace" begin
            ref = zero_state(space)
            discrete = KrylovEvolution(ref, map(_->dt, 0.0:dt:0.2), map(h, 0.0:dt:0.2))
            emulate!(discrete)

            reg = zero_state(space, ComplexLayout())
            prob = SchrodingerProblem(reg, 0.2, h)
            emulate!(prob)
            @test prob.reg ≈ ref atol=1e-4
        end

        @testset "fullspace" begin
            ref = zero_state(5)
            emulate!(KrylovEvolution(ref, map(_->dt, 0.0:dt:0.2), map(h, 0.0:dt:0.2)))

            reg = zero_state(5)
            emulate!(SchrodingerProblem(reg, 0.2, h))
            @test ref ≈ reg atol=1e-4
        end
    end
end

@testset "assertions" begin
    @test_throws ArgumentError SchrodingerProblem(zero_state(10), 0.2, rydberg_h(square_lattice(5, 0.8); Ω=1.0))
end
