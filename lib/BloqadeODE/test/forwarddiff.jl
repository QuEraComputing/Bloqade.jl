using Test
using Random
using ForwardDiff
using BloqadeODE
using YaoSubspaceArrayReg
using FiniteDifferences
using YaoArrayRegister
using BloqadeMIS
using BloqadeWaveforms

function loss(xs::Vector)
    reg = zero_state(Complex{eltype(xs)}, 5)
    tspan = (0, 0.1)
    atoms = [(i,) for i in 1:5]
    Ω = piecewise_constant(clocks = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5], values = xs)
    h = rydberg_h(atoms; C = 2π * 109.16, Ω, ϕ = 0.1)
    prob = SchrodingerProblem(reg, tspan, h; dt = 1e-3)
    emulate!(prob)
    return fidelity(reg, product_state(bit"00001"))
end

Random.seed!(42)
xs = rand(5)

Δ_ad = ForwardDiff.gradient(loss, xs)
Δ_fd, = FiniteDifferences.grad(central_fdm(5, 1), loss, xs)

@test Δ_ad ≈ Δ_fd rtol = 1e-2
