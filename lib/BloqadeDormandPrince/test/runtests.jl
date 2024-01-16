using Test
using BloqadeDormandPrince
using BloqadeExpr
using BloqadeWaveforms
using YaoArrayRegister
using LinearAlgebra




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
    integrate!(BloqadeDPSolver(zero_state(1), (0, 4), h), times) do t, state
        push!(dp_vals, statevec(copy(state)))# register instance returned, get underlying vector out
    end

    @test dp_vals ≈ exact_vals

end

# test the BloqadeDPSolver no longer holds a copy of the register, it just points
# directly to the inputted register which means it mutates the register
@testset "In-Place Register" begin
    atoms = [(i,) for i in 1:5]
    h = rydberg_h(atoms; C = 2π * 109.16, Ω = sin, ϕ = cos)

    # DormandPrince.jl 
    dp_reg = zero_state(5)
    bs = BloqadeDPSolver(dp_reg, (0, 10), h; copy_init=false)
    integrate!(bs, 1.0)

    # strong equality
    @test dp_reg === get_current_state(bs)
    @test dp_reg === register(bs)
    @test emulate!(bs) === bs
end

@testset "Solver Pretty Printing" begin
    atoms = [(i,) for i in 1:10]
    h = rydberg_h(atoms; C = 2π * 109.16, Ω = sin, ϕ = cos)
    solver = BloqadeDPSolver(zero_state(10), 1.5 , h)

    show(stdout, MIME"text/plain"(), solver)
end