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

function to_schema(h::AbstractBlock)
    atoms::Maybe{Vector} = nothing
    # Φ::Maybe{SumOfXPhase} = nothing
    ϕ::Maybe{Waveform} = nothing
    # Ω::Maybe{SumOfX} = nothing
    Ω::Maybe{Waveform} = nothing
    # Δ::Maybe{SumOfN} = nothing
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
        effective_hamiltonian=to_hamiltonian(; ϕ=ϕ, Ω=Ω, Δ=Δ)
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

function to_hamiltonian(; ϕ::Waveform, Ω::Waveform, Δ::Waveform)
    return EffectiveHamiltonian(; rydberg=RydbergHamiltonian(;
        rabi_frequency_amplitude=RydbergRabiFrequencyAmplitude(;
            global_value=RydbergRabiFrequencyAmplitudeGlobal(;
                times=Ω.f.clocks,
                values=Ω.f.values
            )
        ),
        rabi_frequency_phase=RydbergRabiFrequencyPhase(;
            global_value=RydbergRabiFrequencyPhaseGlobal(;
                times=ϕ.f.clocks,
                values=ϕ.f.values
            )
        ),
        detuning=RydbergDetuning(;
            global_value=RydbergDetuningGlobal(;
                times=Δ.f.clocks,
                values=Δ.f.values
            )
        )
    ))
end
