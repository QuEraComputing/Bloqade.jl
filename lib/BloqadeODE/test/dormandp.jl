using Test
using BloqadeODE: BloqadeSolver, DP8Solver, DP5Solver, SolverIterator, integrate!
using BloqadeExpr: rydberg_h, nqudits
using Bloqade
using YaoArrayRegister


@testset "Single End Times" begin
    @testset "DP5" begin
        atoms = [(i,) for i in 1:5]
        h = rydberg_h(atoms; C = 2π * 109.16, Ω = sin, ϕ = cos)

        # DormandPrince.jl 
        dp_reg = zero_state(5)
        bs = BloqadeSolver(dp_reg, 0.0, h; solver_type=DP5Solver)
        integrate!(bs, 1.0)

        # OrdinaryDiffEq 
        tspan = (0, 1.0)
        ode_reg = zero_state(5)
        problem = SchrodingerProblem(ode_reg, tspan, h; algo=DP5())
        emulate!(problem)

        @test dp_reg ≈ ode_reg 
    end
end

h = rydberg_h([(0,0), (0,6)], Ω=t -> sin(t))
reg = zero_state(nqudits(h))


integrate!(BloqadeSolver(reg, 0, h), collect(0:0.01:1.0)) do t, state
    [rydberg_density(state, i) for i in 1:nqudits(h)]
end
