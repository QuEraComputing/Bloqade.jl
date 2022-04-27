# using YaoBlocks
using Yao

function to_json(h::AbstractBlock)
    # 1. check if the input block expression is a summation
    # of RydInteract, SumOfX, SumOfXPhase, SumOfN

    components = Set(map((x) -> typeof(content(x)), h))
    if length(components) > 0 && components <= Set([RydInteract, SumOfX, SumOfXPhase, SumOfN])
        return to_schema(h).to_dict()
    else
        return nothing
    end


    # 2. extract the atom positions, Ω, ϕ, Δ (inside RydbergInteract, )


    # 3. generate corresponding pulses in QuEraSchema
    # 4. call to_dict on QuEraSchema then using JSON3
    #    to convert the dict to json string
    # NOTE: let's ignore local pulse pattern for now
    # since in our representation it's all local pulse
end

function to_schema(h::AbstractBlock, max_slope::Number)
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
        elseif typeof(contents) == SumOfN
            Δ = contents.Δ
        end
    end

    return TaskSpecification(;
        nshots=1,
        lattice=to_lattice(atoms),
        effective_hamiltonian=to_hamiltonian(; ϕ=ϕ, Ω=Ω, Δ=Δ, max_slope=max_slope)
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

# TODO: Check if it's piecewise constant and if so then use the max slope to convert to piecewise linear.
# Given a piecewise constant function with clocks [t1, t2, t3] and values [v1, v2, v3], this creates 
# clocks [t1, t2-((v2-v1)/max_slope), t2, v3-((v3-v2)/max_slope), v3] and values [v1, v1, v2, v2, v3]
function get_piecewise_linear_times_and_clocks(w::Waveform, max_slope::Number)
    isempty(w.f.clocks) && return ([], [])
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

function to_hamiltonian(; ϕ::Waveform, Ω::Waveform, Δ::Waveform, max_slope::Number)
    amp_times, amp_values = get_piecewise_linear_times_and_clocks(Ω, max_slope)
    phase_times, phase_values = get_piecewise_linear_times_and_clocks(ϕ, max_slope)
    detuning_times, detuning_values = get_piecewise_linear_times_and_clocks(Δ, max_slope)

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
