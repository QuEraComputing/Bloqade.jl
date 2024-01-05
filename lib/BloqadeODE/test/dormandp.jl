using Test
using LinearAlgebra

using BloqadeODE: BloqadeSolver, DP8Solver, DP5Solver, SolverIterator, integrate!
using BloqadeExpr: rydberg_h, nqudits
using Bloqade
using YaoArrayRegister


@testset "Single End Time Interface" begin
    atoms = [(i,) for i in 1:5]
    h = rydberg_h(atoms; C = 2π * 109.16, Ω = sin, ϕ = cos)

    # DormandPrince.jl 
    dp_reg = zero_state(5)
    bs = BloqadeSolver(dp_reg, 0, h; solver_type=DP5Solver)
    integrate!(bs, 1.0)

    # OrdinaryDiffEq 
    tspan = (0, 1.0)
    ode_reg = zero_state(5)
    problem = SchrodingerProblem(ode_reg, tspan, h; algo=DP5())
    emulate!(problem)

    @test get_current_state(bs) ≈ ode_reg 
end

@testset "Multiple End Time Interface" begin

    times = [0.1, 0.2, 0.3, 0.4]

    # exact solution
    function evolution_operator(t::Float64)
        ϕ = 2.2 * sin(π * t)^2
        U = zeros(ComplexF64, 2,2)
        U[1,1] =  1 / sqrt(2)
        U[2,1] =  1 / sqrt(2)
        U[2,2] =  1 / sqrt(2)
        U[1,2] = -1 / sqrt(2)
    
        U * diagm(exp.([-im*ϕ, im*ϕ])) * U'
    end
    
    function solution(t)
        U = evolution_operator(t)
        return U * [1.0, 0.0]
    end

    exact_vals = []

    for time in times
        push!(exact_vals, solution(time))
    end

    # equivalent Hamiltonian to exact solution
    atoms = [(0,0)]
    wf = Waveform(t->2.2*2π*sin(2π*t), duration = 1.3);
    h = rydberg_h(atoms; Ω = wf)

    # Bloqade/Dormand Prince multiple end times interface
    dp_vals = []
    integrate!(BloqadeSolver(zero_state(1), 0, h), times) do t, state
        push!(dp_vals, statevec(copy(state)))# register instance returned, get underlying vector out
    end

    @test dp_vals ≈ exact_vals

end

# test the BloqadeSolver no longer holds a copy of the register, it just points
# directly to the inputted register which means it mutates the register
@testset "In-Place Register" begin
    atoms = [(i,) for i in 1:5]
    h = rydberg_h(atoms; C = 2π * 109.16, Ω = sin, ϕ = cos)

    # DormandPrince.jl 
    dp_reg = zero_state(5)
    bs = BloqadeSolver(dp_reg, 0, h; copy_init=false)
    integrate!(bs, 1.0)

    # OrdinaryDiffEq 
    tspan = (0, 1.0)
    ode_reg = zero_state(5)
    problem = SchrodingerProblem(ode_reg, tspan, h)
    emulate!(problem)

    @test dp_reg ≈ ode_reg 
    # strong equality
    @test dp_reg === get_current_state(bs)
end

#=
h = rydberg_h([(0,0), (0,6)], Ω=t -> sin(t))
reg = zero_state(nqudits(h))


integrate!(BloqadeSolver(reg, 0, h), collect(0:0.01:1.0)) do t, state
    [rydberg_density(state, i) for i in 1:nqudits(h)]
end
=#
