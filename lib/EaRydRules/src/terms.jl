function (eq::SchrodingerEquation)(state, p, t::Number)
    dstate = similar(state)
    update_term!(eq.cache.hamiltonian, eq.hamiltonian(t), eq.space)
    mul!(eq.cache.state, eq.cache.hamiltonian, state)
    # @. dstate = -im * eq.cache.state
    update_dstate!(dstate, eq.cache.state, eq.layout)
    return dstate
end

function ChainRulesCore.rrule(eq::SchrodingerEquation, state, p, t::Number)
    dstate = eq(state, p, t)
    function pullback(Δ)
        # s̄ = iH * Δ
        # ̅Ω = ℜ[Δ' * -i (∑ σₓ/2) * state]
        # ̅Δ = ℜ[Δ' * -i (∑ n) * state]
        # c̄ᵢⱼ = ℜ[Δ' * -i (nᵢnⱼ) * state]
        s̄ = -Δ
        update_term!(eq.cache.hamiltonian, eq.hamiltonian(t), eq.space)
        mul!(eq.cache.state, eq.cache.hamiltonian, s̄)
        update_dstate!(s̄, eq.cache.state, eq.layout)
        _iΔconj = -im*Δ'
        gterms = map(term->back_term(_iΔconj, term, state, eq.space), terms)
        ē = Tangent{typeof(eq)}(; layout=NoTangent(),
            hamiltonian=Tangent{typeof(eq.hamiltonian)}(; terms= (gΩ, gΔ, gcs)),
            space=NoTangent(), cache=NoTangent())
        return ē, s̄, NoTangent()
    end
    return dstate, pullback
end

function back_term(_iΔconj, term::ZTerm, state, space)
    gΔs = zero(term.Δs)
    for k = 1:length(state)
        b = space[k]
        for i=1:term.nsites
            gΔs[i] += _iΔconj[k] * state[k] * (b >> (i-1) & 1)
        end
    end
    return gΔs
end

function back_term(_iΔconj, term::XTerm, state, space)
    gϕs = zero(term.ϕs)
    gΩs = zero(term.Ωs)
    for k = 1:length(state)
        b = space[k]
        for i=1:term.nsites
            gΔs[i] += _iΔconj[k] * state[k] * (b >> (i-1) & 1)
        end
    end
    return gΔs
end

