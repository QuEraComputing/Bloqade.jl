function rydberg_h(atom_positions; C::Real, Ω, ϕ, Δ)
    positions = map(atom_positions) do position
        (position..., )
    end

    nsites = length(positions)
    term = RydInteract(positions, C)

    Ω = div_by_two(Ω)

    if !isnothing(Ω) && !isnothing(ϕ)
        term += SumOfXPhase(nsites, Ω, ϕ)
    elseif !isnothing(Ω) && isnothing(ϕ)
        term += SumOfX(nsites, Ω)
    end

    if !isnothing(Δ)
        term -= SumOfN(nsites, Δ)
    end

    return term
end

function div_by_two(Ω)
    isnothing(Ω) && return
    if !is_time_function(Ω)
        return Ω ./ 2
    end

    return if Ω isa Vector
        map(Ω) do Ω_i
            t->Ω_i(t) / 2
        end
    else
        t->Ω(t)/2
    end
end
