using Test
using ForwardDiff
using BloqadeExpr
using YaoArrayRegister
using BloqadeWaveforms
using BloqadeKrylov
using BloqadeLattices

function test(x;tol=1e-7)
    T = 0.1
    dT = T/10
    reg = zero_state(Complex{typeof(x)}, 5)
    atoms = [(6.7*i,) for i in 1:5]
    h = rydberg_h(atoms; Ω = x, Δ = t->-5*cos(π*t/T))
    prob = KrylovEvolution(reg, 0:dT:T, h;tol=tol)
    emulate!(prob)
    return abs2.(statevec(reg))
end

test(10.0;tol=1e-9)


# ForwardDiff.derivative(2.0) do x
#     reg = zero_state(Complex{typeof(x)}, 5)
#     atoms = [(i,) for i in 1:5]
#     h = rydberg_h(atoms; Ω = sin, Δ = x)
#     prob = KrylovEvolution(reg, 0:1e-3:0.1, h)
#     emulate!(prob)
#     return abs2(statevec(reg)[1])
# end

# using ExponentialUtilities

# ForwardDiff.derivative(2.0) do x
#     st = rand(2)
#     st = BloqadeKrylov.expmv(2.0, fill(x, (2, 2)), st)
#     return sum(st)
# end
