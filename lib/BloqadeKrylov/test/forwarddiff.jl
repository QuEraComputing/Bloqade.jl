using Test
using ForwardDiff
using BloqadeExpr
using YaoArrayRegister
using BloqadeKrylov


ForwardDiff.derivative(2.0) do x
    T = 0.1
    dT = T / 10
    reg = zero_state(Complex{typeof(x)}, 5)
    atoms = [(6.7 * i,) for i in 1:5]
    h = rydberg_h(atoms; Ω = x, Δ = t -> -5 * cos(π * t / T))
    prob = KrylovEvolution(reg, 0:dT:T, h)
    emulate!(prob)
    return abs2(statevec(reg)[1])
end

using ExponentialUtilities

@testset "Gradients of atom positions" begin
    nsites = 5
    scale = 5.0
    function loss_atoms(ps)
        atoms = [[ps[i, 1], ps[i, 2]] for i = 1:nsites]
        reg = zero_state(Complex{typeof(ps[1, 1])}, nsites)
        rh = rydberg_h(atoms; Ω = 1)
        prob = KrylovEvolution(reg, 0:1e-2:1, rh)
        emulate!(prob)
        return abs2(statevec(reg)[1])
    end
    ps = rand(nsites, 2) * scale
    ForwardDiff.gradient(loss_atoms, ps)
end

@testset "ForwardDiff broken" begin
    ForwardDiff.derivative(2.0) do x
        reg = zero_state(Complex{typeof(x)}, 5)
        atoms = [(6*i,) for i in 1:5]
        h = rydberg_h(atoms; Ω = sin, Δ = x)
        prob = KrylovEvolution(reg, 0:1e-3:0.02, h)
        emulate!(prob)
        return abs2(statevec(reg)[1])
    end
    
    # # failing function, doesn't seem to 
    # @test isnan(ForwardDiff.derivative(2.0) do x
    #     h = fill(zero(typeof(x)), (2, 2))
    #     F = ExponentialUtilities.exponential!(im * h, ExpMethodGeneric())
    #     return sum(F)
    # end)
end
