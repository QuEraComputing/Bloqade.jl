using Yao
using BloqadeODE
using BloqadeExpr
using BloqadeMIS
using BloqadeKrylov
using Test

@testset "ODE problem interface" begin
    reg = zero_state(5)
    tspan = (0, 1e-4)
    atoms = [(i,) for i in 1:5]
    h = rydberg_h(atoms; C = 2π * 109.16, Ω = sin, ϕ = cos)
    prob = SchrodingerProblem(reg, tspan, h; dt = 1e-5, progress = true, save_start = false)

    println()
    show(stdout, MIME"text/plain"(), prob)
    println()
    show(stdout, MIME"text/plain"(), prob.f.f)

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
    atoms = [(i,) for i in 1:5]
    h = rydberg_h(atoms; C = 2π * 109.16, Ω = sin, ϕ = cos)
    prob = SchrodingerProblem(reg, tspan, h; dt = 1e-5, progress = true, save_start = false)
    integrator = init(prob, Vern8())
    for (u, t) in TimeChoiceIterator(integrator, [0.05, 0.1])
        @test t in [0.05, 0.1]
    end
end

@testset "ODE evolution results" begin
    atoms = [(i,) for i in 1:5]
    space = blockade_subspace(atoms, 1.5)

    @testset "h=$name" for (name, h) in [
        "x+z" => SumOfX(5, 1.0) + SumOfN(5, sin),
        "rydberg" => rydberg_h(atoms; Δ = sin, Ω = cos, C = 2π * 109),
    ]
        dt = 1e-5
        @testset "subspace" begin
            ref = zero_state(space)
            discrete = KrylovEvolution(ref, collect(0.0:dt:0.2), h)
            BloqadeExpr.emulate!(discrete)

            reg = zero_state(space)
            prob = SchrodingerProblem(reg, 0.2, h)
            BloqadeExpr.emulate!(prob)
            @test prob.reg ≈ ref atol = 1e-4
        end

        @testset "fullspace" begin
            ref = zero_state(5)
            BloqadeExpr.emulate!(KrylovEvolution(ref, collect(0.0:dt:0.2), h))

            reg = zero_state(5)
            BloqadeExpr.emulate!(SchrodingerProblem(reg, 0.2, h))
            @test ref ≈ reg atol = 1e-4
        end
    end
end

@testset "assertions" begin
    atoms = [(i,) for i in 1:5]
    @test_throws ArgumentError SchrodingerProblem(zero_state(10), 0.2, rydberg_h(atoms; Ω = 1.0))
end

@testset "3-level W-state ODE simulation" begin
    nsites = 2
    V = 1e5
    atoms = [(0.0, 0.0), (0.0, (862690*2pi/V)^(1/6))]
    reg = product_state(dit"11;3")
    h = rydberg_h_3(atoms; Ω_r = 1)

    goal = exp(-im*pi/sqrt(2)*Matrix(h)) * state(reg)
    prob = SchrodingerProblem(reg, (0, pi/sqrt(2)), h)
    emulate!(prob)
    @test isapprox(state(reg), goal; atol = 1e-3)
end
