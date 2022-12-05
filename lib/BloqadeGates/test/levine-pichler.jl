using BloqadeGates
using BloqadeODE
using YaoBlocks
using Test

nsites = 2
V = 1e3
Ω_r = 1.0
τ = 4.29268/Ω_r
Δ_r = 0.377371*Ω_r
ϕ_r = 3.90242  # sign flipped
atoms = [(0.0, 0.0), (0.0, (862690*2pi/V)^(1/6))]

reg = zero_state(nsites; nlevel = 3)
p1 = global_single_qubit_gate(atoms, H; backend = SchrodingerProblem)
reg |> p1
reg |> state

p2 = global_pulse(atoms, Ω_r, 0.0, Δ_r, τ; pulse_type = :rydberg, backend = SchrodingerProblem)
reg |> p2
reg |> state

p3 = global_pulse(atoms, Ω_r, ϕ_r, Δ_r, τ; pulse_type = :rydberg, backend = SchrodingerProblem)
reg |> p3
reg |> state

st = state(reg)[[1,2,4,5]]
st /= st[1]

mat(p2)