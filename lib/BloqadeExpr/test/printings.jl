using Test
using BloqadeExpr

atoms = [(1, 2), (2, 3)]

for h in [
    rydberg_h(atoms; Ω = 1.0, Δ = sin),
    rydberg_h(atoms; Ω = 1.0, Δ = [sin, sin]),
    rydberg_h(atoms; Ω = sin, Δ = 1.0),
    rydberg_h(atoms; Ω = [sin, sin], Δ = 1.0),
    rydberg_h(atoms; Ω = 2.0, Δ = sin),
    rydberg_h(atoms; Ω = 2.0, Δ = [sin, sin]),
    rydberg_h(atoms; Ω = sin, Δ = 2.0),
    rydberg_h(atoms; Ω = [sin, sin], Δ = 2.0),
    rydberg_h(atoms; Ω = [1.0, 2.0], Δ = sin),
    rydberg_h(atoms; Ω = [1.0, 2.0], Δ = [sin, sin]),
    rydberg_h(atoms; Ω = sin, Δ = [1.0, 2.0]),
    rydberg_h(atoms; Ω = [sin, sin], Δ = [1.0, 2.0]),
    rydberg_h(atoms; Ω = 1.0, ϕ = sin),
    rydberg_h(atoms; Ω = 1.0, ϕ = [sin, sin]),
    rydberg_h(atoms; Ω = sin, ϕ = 1.0),
    rydberg_h(atoms; Ω = [sin, sin], ϕ = 1.0),
    rydberg_h(atoms; Ω = 2.0, ϕ = sin),
    rydberg_h(atoms; Ω = 2.0, ϕ = [sin, sin]),
    rydberg_h(atoms; Ω = sin, ϕ = 2.0),
    rydberg_h(atoms; Ω = [sin, sin], ϕ = 2.0),
    rydberg_h(atoms; Ω = [1.0, 2.0], ϕ = sin),
    rydberg_h(atoms; Ω = [1.0, 2.0], ϕ = [sin, sin]),
    rydberg_h(atoms; Ω = sin, ϕ = [1.0, 2.0]),
    rydberg_h(atoms; Ω = [sin, sin], ϕ = [1.0, 2.0]),
]
    show(stdout, MIME"text/plain"(), h)
    println(stdout, "----------")
    println(stdout, "----------")
    show(stdout, MIME"text/latex"(), h)

    show(stdout, MIME"text/plain"(), Hamiltonian(ComplexF64, h))
end

params = [nothing, 1.0, [1.0, 2.0], sin, [sin, cos]]
for Ω_hf in params, ϕ_hf in params
    h = rydberg_h_3(atoms; Ω_hf, Δ_hf = 1.0, ϕ_hf, Ω_r = 1.0, Δ_r = 1.0, ϕ_r = 1.0)
    show(stdout, MIME"text/plain"(), h)
    println(stdout, "----------")
    println(stdout, "----------")
    show(stdout, MIME"text/latex"(), h)
    show(stdout, MIME"text/plain"(), Hamiltonian(ComplexF64, h))
end
for Δ_hf in params
    h = rydberg_h_3(atoms; Ω_hf = 1.0, Δ_hf, ϕ_hf = 1.0, Ω_r = 1.0, Δ_r = 1.0, ϕ_r = 1.0)
    show(stdout, MIME"text/plain"(), h)
    println(stdout, "----------")
    println(stdout, "----------")
    show(stdout, MIME"text/latex"(), h)
    show(stdout, MIME"text/plain"(), Hamiltonian(ComplexF64, h))
end
for Ω_r in params, ϕ_r in params
    h = rydberg_h_3(atoms; Ω_hf = 1.0, Δ_hf = 1.0, ϕ_hf = 1.0, Ω_r, Δ_r = 1.0, ϕ_r)
    show(stdout, MIME"text/plain"(), h)
    println(stdout, "----------")
    println(stdout, "----------")
    show(stdout, MIME"text/latex"(), h)
    show(stdout, MIME"text/plain"(), Hamiltonian(ComplexF64, h))
end
for Δ_r in params
    h = rydberg_h_3(atoms; Ω_hf = 1.0, Δ_hf = 1.0, ϕ_hf = 1.0, Ω_r = 1.0, Δ_r, ϕ_r = 1.0)
    show(stdout, MIME"text/plain"(), h)
    println(stdout, "----------")
    println(stdout, "----------")
    show(stdout, MIME"text/latex"(), h)
    show(stdout, MIME"text/plain"(), Hamiltonian(ComplexF64, h))
end

for op in (SumOfZ, SumOfZ_01, SumOfZ_1r), p in params
    h = op(2, isnothing(p) ? 0 : p)
    println(stdout, "\n")
    show(stdout, MIME"text/plain"(), h)
    println("")
    show(stdout, MIME"text/latex"(), h)
    println(stdout, "\n----------")
    println(stdout, "----------")
    show(stdout, MIME"text/plain"(), Hamiltonian(ComplexF64, h))
end