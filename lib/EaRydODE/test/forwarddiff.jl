using Test
using ForwardDiff
using EaRydODE
using YaoSubspaceArrayReg
using YaoArrayRegister

ForwardDiff.derivative(2.0) do x
    reg = zero_state(Complex{typeof(x)}, 5)
    tspan = (0, 1e-4)
    atoms = [(i, ) for i in 1:5]
    h = rydberg_h(atoms; C=2π * 109.16, Ω=sin, ϕ=x)
    prob = SchrodingerProblem(reg, tspan, h; dt=1e-5, progress=true, save_start=false)
    emulate!(prob)
    sum(abs2.(statevec(prob.reg)))
end
