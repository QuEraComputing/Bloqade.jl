# using YaoBlocks
using Yao
using Configurations

function to_json(h::AbstractBlock, params::SchemaConversionParams)
    # 1. check if the input block expression is a summation
    # of RydInteract, SumOfX, SumOfXPhase, SumOfN

    # components = Set(map((x) -> typeof(content(x)), h))
    # if length(components) > 0 && components <= Set([RydInteract, SumOfX, SumOfXPhase, SumOfN])
    #     return to_schema(h).to_dict()
    # else
    #     return nothing
    # end
    return Configurations.to_dict(to_schema(h;
        rabi_frequency_amplitude_max_slope=params.rabi_frequency_amplitude_max_slope,
        rabi_frequency_phase_max_slope=params.rabi_frequency_phase_max_slope,
        rabi_detuning_max_slope=params.rabi_detuning_max_slope,
        n_shots=params.n_shots
    ))


    # 2. extract the atom positions, Ω, ϕ, Δ (inside RydbergInteract, )


    # 3. generate corresponding pulses in QuEraSchema
    # 4. call to_dict on QuEraSchema then using JSON3
    #    to convert the dict to json string
    # NOTE: let's ignore local pulse pattern for now
    # since in our representation it's all local pulse
end

function to_schema(h::AbstractBlock; rabi_frequency_amplitude_max_slope::Number,
    rabi_frequency_phase_max_slope::Number, rabi_detuning_max_slope::Number, n_shots::Number
)
    atoms::Maybe{Vector} = nothing
    ϕ::Maybe{Waveform} = nothing
    Ω::Maybe{Waveform} = nothing
    Δ::Maybe{Waveform} = nothing

    for component in h
        contents = content(component)
        if typeof(contents) == RydInteract
            atoms = contents.atoms
        elseif typeof(contents) == SumOfXPhase
            ϕ = contents.ϕ
            Ω = contents.Ω.Ω
        elseif typeof(contents) == SumOfX
            Ω = contents.Ω.Ω
        elseif typeof(contents) == SumOfN
            Δ = contents.Δ
        end
    end

    return TaskSpecification(;
        nshots=n_shots,
        lattice=to_lattice(atoms),
        effective_hamiltonian=to_hamiltonian(; ϕ=ϕ, Ω=Ω, Δ=Δ,
            rabi_frequency_amplitude_max_slope=rabi_frequency_amplitude_max_slope,
            rabi_frequency_phase_max_slope=rabi_frequency_phase_max_slope,
            rabi_detuning_max_slope=rabi_detuning_max_slope
        )
    )

end

function to_lattice(atoms::Vector)
    coords = Vector{Tuple{Float64,Float64}}()
    for atom in atoms
        coord = atom
        if length(atom) == 1
            coord = Tuple([atom[1], 0])
        end
        push!(coords, coord)
    end
    return Lattice(; sites=coords, filling=vec(ones(length(coords), 1)))
end

# Check if the waveform is piecewise constant and if so then use the max slope to convert to piecewise linear.
# Given a piecewise constant function with clocks [t1, t2, t3] and values [v1, v2, v3], this creates 
# clocks [t1, t2-((v2-v1)/max_slope), t2, v3-((v3-v2)/max_slope), v3] and values [v1, v1, v2, v2, v3]
# If the waveform is empty or nothing, returns clocks [0] and values [0]
function get_piecewise_linear_times_and_clocks(w::Maybe{Waveform}, max_slope::Number)
    # if nothing or empty
    (isnothing(w) || isempty(w.f.clocks)) && return ([0], [0])

    # if not piecewise constant, return clocks and values directly
    !isa(w.f, BloqadeWaveforms.PiecewiseConstant) && return (w.f.clocks, w.f.values)

    clocks = []
    values = []
    for i in 1:(length(w.f.clocks)-1)
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

function to_hamiltonian(; ϕ::Maybe{Waveform}, Ω::Maybe{Waveform}, Δ::Maybe{Waveform}, rabi_frequency_amplitude_max_slope::Number,
    rabi_frequency_phase_max_slope::Number, rabi_detuning_max_slope::Number
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
