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


@testset "ForwardDiff broken" begin
    # # failing function, doesn't seem to 
    # ForwardDiff.derivative(2.0) do x
    #     reg = zero_state(Complex{typeof(x)}, 5)
    #     atoms = [(6*i,) for i in 1:5]
    #     h = rydberg_h(atoms; Ω = sin ,Δ = x)
    #     prob = KrylovEvolution(reg, 0:1e-3:0.02, h)
    #     emulate!(prob)
    #     return abs2(statevec(reg)[1])
    # end

    @test isnan(ForwardDiff.derivative(2.0) do x
        h = fill(zero(typeof(x)), (2, 2))
        F = ExponentialUtilities.exponential!(im * h, ExpMethodGeneric())
        return sum(F)
    end)

end
