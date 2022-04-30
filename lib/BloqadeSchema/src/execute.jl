using Yao
using Configurations
using JSON3

function to_json(h::AbstractBlock, params::SchemaConversionParams)
    return JSON3.write(BloqadeSchema.to_dict(h, params))
end

function to_dict(h::AbstractBlock, params::SchemaConversionParams)
    return Configurations.to_dict(to_schema(h;
        rabi_frequency_amplitude_max_slope=params.rabi_frequency_amplitude_max_slope,
        rabi_frequency_phase_max_slope=params.rabi_frequency_phase_max_slope,
        rabi_detuning_max_slope=params.rabi_detuning_max_slope,
        n_shots=params.n_shots
    ))
end

function assert_hamiltonian_schema(h::AbstractBlock)
    h isa RydInteract || return
    h isa Add || error("expect a Rydberg hamiltonian")
    isempty(blockfilter(x->isa(x, RydInteract), h)) && error("expect RydInteract term")

    for each in subblocks(h)
        if each isa Scale
            content(each) isa Union{RydInteract, SumOfX, SumOfXPhase} && factor(each) == 1 ||
                error("only hamiltonian created by rydberg_h is supported")
            content(each) isa SumOfN && factor(each) == -1 || error("expect the prefactor of SumOfN to be -1")
        end
    end
    return
end

function to_schema(h::AbstractBlock; rabi_frequency_amplitude_max_slope::Real,
    rabi_frequency_phase_max_slope::Real, rabi_detuning_max_slope::Real, n_shots::Real
)
    assert_hamiltonian_schema(h)
    atoms = nothing
    ϕ = nothing
    Ω = nothing
    Δ = nothing

    for each in subblocks(h)
        block = content(each)
        if block isa RydInteract
            atoms = block.atoms
        elseif block isa SumOfXPhase
            ϕ = block.ϕ
            Ω = block.Ω.Ω
        elseif block isa SumOfX
            Ω = block.Ω.Ω
        elseif block isa SumOfN
            Δ = block.Δ
        end
    end

    return TaskSpecification(;
        nshots=n_shots,
        lattice=to_lattice(atoms),
        effective_hamiltonian=to_hamiltonian(Ω, ϕ, Δ,
            rabi_frequency_amplitude_max_slope,
            rabi_frequency_phase_max_slope,
            rabi_detuning_max_slope
        )
    )

end

function to_lattice(atoms::Vector)
    coords = map(atoms) do coord
        length(coord) == 1 && return (coord[1], 0)
        coord
    end
    return Lattice(; sites=coords, filling=vec(ones(length(coords), 1)))
end

# Check if the waveform is piecewise constant and if so then use the max slope to convert to piecewise linear.
# Given a piecewise constant function with clocks [t1, t2, t3] and values [v1, v2, v3], this creates 
# clocks [t1, t2-((v2-v1)/max_slope), t2, v3-((v3-v2)/max_slope), v3] and values [v1, v1, v2, v2, v3]
# If the waveform is empty or nothing, returns clocks [0] and values [0]
function get_piecewise_linear_times_and_clocks(w::Maybe{Waveform}, max_slope::Real)
    # if nothing or empty
    (isnothing(w) || isempty(w.f.clocks)) && return ([0], [0])

    # if not piecewise constant, return clocks and values directly
    !isa(w.f, BloqadeWaveforms.PiecewiseConstant) && return (w.f.clocks, w.f.values)

    clocks = Real[]
    values = Real[]
    for i in 1:(length(w.f.values)-1)
        rise = abs(w.f.values[i+1] - w.f.values[i])
        run_t = rise / max_slope

        append!(clocks, w.f.clocks[i])
        append!(values, w.f.values[i])

        run_t != 0 && append!(clocks, w.f.clocks[i+1] - run_t)
        run_t != 0 && append!(values, w.f.values[i])
    end

    append!(clocks, last(w.f.clocks))
    append!(values, last(w.f.values))
    return (clocks, values)

end

function to_hamiltonian(Ω::Maybe{Waveform}, ϕ::Maybe{Waveform}, Δ::Maybe{Waveform},
    rabi_frequency_amplitude_max_slope::Real, rabi_frequency_phase_max_slope::Real, rabi_detuning_max_slope::Real
)
    amp_times, amp_values = get_piecewise_linear_times_and_clocks(Ω, rabi_frequency_amplitude_max_slope)
    phase_times, phase_values = get_piecewise_linear_times_and_clocks(ϕ, rabi_frequency_phase_max_slope)
    detuning_times, detuning_values = get_piecewise_linear_times_and_clocks(Δ, rabi_detuning_max_slope)

    return EffectiveHamiltonian(; rydberg=RydbergHamiltonian(;
        rabi_frequency_amplitude=RydbergRabiFrequencyAmplitude(;
            global_value=RydbergRabiFrequencyAmplitudeGlobal(;
                times=amp_times,
                values=amp_values
            )
        ),
        rabi_frequency_phase=RydbergRabiFrequencyPhase(;
            global_value=RydbergRabiFrequencyPhaseGlobal(;
                times=phase_times,
                values=phase_values
            )
        ),
        detuning=RydbergDetuning(;
            global_value=RydbergDetuningGlobal(;
                times=detuning_times,
                values=detuning_values
            )
        )
    ))
end
