using Test
using ForwardDiff
using BloqadeExpr
using YaoArrayRegister
using BloqadeWaveforms
using BloqadeKrylov
using BloqadeLattices


ForwardDiff.derivative(2.0) do x
    reg = zero_state(Complex{typeof(x)}, 5)
    atoms = [(i, ) for i in 1:5]
    h = rydberg_h(atoms; Ω=sin, Δ=x)
    prob = KrylovEvolution(reg, 0:1e-3:0.1, h)
    emulate!(prob)
    abs2(statevec(reg)[1])
end

using ExponentialUtilities

ForwardDiff.derivative(2.0) do x
    st = rand(2)
    st = BloqadeKrylov.expmv(2.0, fill(x, (2, 2)), st)
    sum(st)
end
