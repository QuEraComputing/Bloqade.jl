


function mult_by_two(Ω)
    isnothing(Ω) && return
    if BloqadeExpr.is_const_param(Ω)
        return Ω .* 2
    end

    return if Ω isa Vector
        map(Ω) do Ω_i
            return Ω_i.f
        end
    else
        Ω.f
    end

end


function get_rydberg_params(h::BloqadeExpr.RydbergHamiltonian)
    # extracts parameters from RydbergHamiltonian
    ϕ = nothing
    Ω = nothing
    Δ = nothing

    if h.rabi_term isa BloqadeExpr.SumOfX
        Ω = mult_by_two(h.rabi_term.Ω)
    elseif h.rabi_term isa BloqadeExpr.SumOfXPhase
        Ω = mult_by_two(h.rabi_term.Ω)
        ϕ = h.rabi_term.ϕ
    end

    if h.detuning_term isa BloqadeExpr.SumOfN
        Δ = h.detuning_term.Δ
    end

    return (h.rydberg_term.atoms,ϕ,Ω,Δ)
end


schema_parse_ϕ(ϕ::Waveform{BloqadeWaveforms.PiecewiseLinear{T,I},T}) where {T<:Real,I} = ϕ
schema_parse_Ω(Ω::Waveform{BloqadeWaveforms.PiecewiseLinear{T,I},T}) where {T<:Real,I} = Ω

schema_parse_ϕ(_) = error("Cannot convert Hamiltonian to schema, ϕ(t) must be piecewise linear Waveform.")
schema_parse_Ω(_) = error("Cannot convert Hamiltonian to schema, Ω(t) must be piecewise linear Waveform.")

function schema_parse_Δ(Δ)
    Δ isa PiecewiseLinearWaveform  && return (Δ,nothing,1.0)
    
    if Δ isa Vector

        map(Δ) do ele
            ele isa PiecewiseLinearWaveform || error("Cannot convert Hamiltonian to schema, Δ(t) must be a single or collection of piecewise linear Waveform(s).")
        end

        clocks = Δ[1].f.clocks
        map(Δ) do ele
            clocks == ele.f.clocks || error("Cannot convert Hamiltonian to schema, all piecewise linear Waveform(s) must have same clocks.")
        end

        Δti = zeros(length(clocks),length(Δ))
    
        for (j,δ) in enumerate(Δ)
            Δti[:,j] .= δ.f.values # copy values to Δti
        end
        
        ((δ_values,Δi),(Δ_values,_)) = find_local_masks(Δti;name=:Δ,assert_truncation=true)
        
        δ_values = round.(δ_values;sigdigits=14)
        Δ_values = round.(Δ_values;sigdigits=14)
        Δi = round.(Δi;sigdigits=14)
        
        
        # get waveforms
        Δ = piecewise_linear(;
            clocks = clocks,
            values = Δ_values
        )

        δ = piecewise_linear(;
            clocks = clocks,
            values = δ_values
        )

        return (Δ,δ,Δi)
    else
        error("Cannot convert Hamiltonian to schema, Δ(t) must be a single or collection of piecewise linear Waveform(s).")
    end
end

# pass fail function. 
# certain criteria is required to pass
# 1. all waveforms must be PWL
# 2. Local detuning must be PWL with all same clocks and can be decomposed 
#    with SVD such that the waveform is Δ(i,t) = Δ(t) + Δ_i * δ(t)
function schema_parse_rydberg_fields(H::BloqadeExpr.RydbergHamiltonian)
    atoms,ϕ,Ω,Δ = get_rydberg_params(H)
    ϕ = schema_parse_ϕ(ϕ)
    Ω = schema_parse_Ω(Ω)
    Δ,δ,Δi = schema_parse_Δ(Δ)
    return (atoms,ϕ,Ω,Δ,δ,Δi)
end
    


