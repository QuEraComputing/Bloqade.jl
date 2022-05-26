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
