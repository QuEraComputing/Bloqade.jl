using EaRydRules
using Test
using OrdinaryDiffEq
using Zygote

@testset "rules" begin
    C, Ω, ϕ, Δ = 200.0, 3.0, 2.0, 1.4
    function f(h::Hamiltonian)
        loss = 0.0
        for term in terms
            if term isa RydInteract
                loss += term.C * 2
                for atom in term.atoms
                    loss += atom[1] * 2
                    loss += atom[2] * 2
                end
            elseif term isa XTerm
                for Ω in term.Ωs
                    loss += Ω * 3
                end
                for ϕ in term.ϕs
                    loss += ϕ * 3
                end
            elseif term isa ZTerm
                for Δ in term.Δs
                    loss += Δ * 4
                end
            else
                error("")
            end
        end
    end
    h, back = Zygote.pullback(rydberg_h, atoms, C, Ω, ϕ, Δ)
    #gatoms, gC, gΩ, gϕ, gΔ = back(Tangent{})
    #RydInteract(atoms, C) + XTerm(length(atoms), Ω, ϕ) - NTerm(length(atoms), Δ)
end

@testset "gradients" begin
    u0 = [1.0]
    f(x, p, t) = [sin((@show x[1]) * p[1])/t]
    tspan = [1, 2.0]
    p = [3.0]
    prob = ODEProblem(f, u0, tspan, p)
    integrator = init(prob, Vern8(), save_everystep=true)
    step!(integrator, 0.1, true)
    @test integrator.t == 1.1
    for i in integrator
    end
    @show integrator
    #treeverse(f, gf, state)
end