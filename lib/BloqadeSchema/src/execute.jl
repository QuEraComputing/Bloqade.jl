# """
#     execute(j::String)

# Executes a task given as a JSON string in the task specification API format, and returns a JSON string of the result
# """
# function execute(j::String)
#     task = Configurations.from_dict(BloqadeSchema.TaskSpecification, JSON.parse(j))
#     h = from_json(j)
#     return JSON.json(Configurations.to_dict(execute(h, length(content(h[1]).atoms), 3e-6, task.nshots)))
# end

# function execute(h::Add, n_atoms::Int, total_time::Float64, nshots::Int)
#     # Always start off with the zero state
#     reg = zero_state(n_atoms)
#     problem = SchrodingerProblem(reg, total_time, h)
#     emulate!(problem)
#     bitstrings = reg |> measure(; nshots=nshots)
#     return to_task_output(bitstrings)
# end

function to_task_output(bitstrings::Vector{<:BitBasis.BitStr64})
    shot_outputs = map(bitstrings) do bs
        # Assume perfect loading/sorting, so the initial loading is full
        pre_sequence = ones(length(bitstrings[1]))

        post_sequence = []
        for i in 1:length(bitstrings[1])
            append!(post_sequence, Int32(readbit(bs, i)))
        end

        return ShotOutput(;
            shot_status_code=200,
            pre_sequence=pre_sequence,
            post_sequence=post_sequence
        )
    end

    return TaskOutput(;
        task_status_code=200,
        shot_outputs=shot_outputs
    )
end

function from_json(j::String)
    t = Configurations.from_dict(BloqadeSchema.TaskSpecification, JSON.parse(j))
    return from_schema(t)
end

function from_schema(t::TaskSpecification)
    atoms = [t.lattice.sites[i] for i in 1:length(t.lattice.sites) if t.lattice.filling[i] == 1]

    rabi_freq_amp = t.effective_hamiltonian.rydberg.rabi_frequency_amplitude.global_value
    rabi_freq_phase = t.effective_hamiltonian.rydberg.rabi_frequency_phase.global_value
    detuning_global = t.effective_hamiltonian.rydberg.detuning.global_value

    Ω = BloqadeWaveforms.piecewise_linear(; clocks=rabi_freq_amp.times, values=rabi_freq_amp.values)
    ϕ = BloqadeWaveforms.piecewise_linear(; clocks=rabi_freq_phase.times, values=rabi_freq_phase.values)
    Δ = BloqadeWaveforms.piecewise_linear(; clocks=detuning_global.times, values=detuning_global.values)

    return BloqadeExpr.rydberg_h(atoms; Δ=Δ, Ω=Ω, ϕ=ϕ)
end

"""
    to_json(h::AbstractBlock; kw...)

Convert a hamiltonian to JSON task specification.

# Examples

```jldoctest
julia> Ω = BloqadeWaveforms.piecewise_constant(; clocks=[0, 2, 4, 6, 7], values=[5, 3, 4, 6]);

julia> Δ = BloqadeWaveforms.piecewise_linear(; clocks=[0.0, 0.6, 2.1, 2.2], values=[-10.1, -10.1, 10.1, 10.1]);

julia> ϕ = BloqadeWaveforms.piecewise_linear(; clocks=[0, 5], values=[33, 0]);

julia> atoms = [(0, 0), (1, 3), (4, 2), (6, 3), (0, 5), (2, 5)];

julia> block = BloqadeExpr.rydberg_h(atoms; Δ=Δ, Ω=Ω, ϕ=ϕ);

julia> BloqadeSchema.to_json(block; n_shots=10)
"{\"nshots\":10,\"lattice\":{\"sites\":[[0.0,0.0],[1.0,3.0],[4.0,2.0],[6.0,3.0],[0.0,5.0],[2.0,5.0]],\"filling\":[1,1,1,1,1,1]},\"effective_hamiltonian\":{\"rydberg\":{\"rabi_frequency_amplitude\":{\"global\":{\"times\":[0.0,-18.0,2.0,-6.0,4.0,-14.0,7.0],\"values\":[5.0,5.0,3.0,3.0,4.0,4.0,6.0]}},\"rabi_frequency_phase\":{\"global\":{\"times\":[0.0,5.0],\"values\":[33.0,0.0]}},\"detuning\":{\"global\":{\"times\":[0.0,0.6,2.1,2.2],\"values\":[-10.1,-10.1,10.1,10.1]}}}}}"
```
"""
to_json(h::AbstractBlock; kw...) = to_json(h, SchemaConversionParams(;kw...))
to_dict(h::AbstractBlock; kw...) = to_dict(h, SchemaConversionParams(;kw...))
to_schema(h::AbstractBlock; kw...) = to_schema(h, SchemaConversionParams(;kw...))

function to_json(h::AbstractBlock, params::SchemaConversionParams)
    return JSON.json(BloqadeSchema.to_dict(h, params))
end

function to_dict(h::AbstractBlock, params::SchemaConversionParams)
    return Configurations.to_dict(to_schema(h, params))
end

function assert_hamiltonian_schema(h::AbstractBlock)
    h isa RydInteract || return
    h isa Add || error("expect a Rydberg hamiltonian")
    isempty(blockfilter(x -> isa(x, RydInteract), h)) && error("expect RydInteract term")

    for each in subblocks(h)
        if each isa Scale
            content(each) isa Union{RydInteract,SumOfX,SumOfXPhase} && factor(each) == 1 ||
                error("only hamiltonian created by rydberg_h is supported")
            content(each) isa SumOfN && factor(each) == -1 || error("expect the prefactor of SumOfN to be -1")
        end
    end
    return
end

function to_schema(h::AbstractBlock, params::SchemaConversionParams)
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
            Ω = block.Ω.f
        elseif block isa SumOfX
            Ω = block.Ω.f
        elseif block isa SumOfN
            Δ = block.Δ
        end
    end

    return TaskSpecification(;
        nshots=params.n_shots,
        lattice=to_lattice(atoms),
        effective_hamiltonian=to_hamiltonian(Ω, ϕ, Δ,
            params.rabi_frequency_amplitude_max_slope,
            params.rabi_frequency_phase_max_slope,
            params.rabi_detuning_max_slope
        )
    )
end

function to_lattice(atoms::Vector)
    coords = map(atoms) do coord
        length(coord) == 1 && return (coord[1], 0)
        return coord
    end
    return Lattice(; sites = coords, filling = vec(ones(length(coords), 1)))
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

function to_hamiltonian(
    Ω::Maybe{Waveform},
    ϕ::Maybe{Waveform},
    Δ::Maybe{Waveform},
    rabi_frequency_amplitude_max_slope::Real,
    rabi_frequency_phase_max_slope::Real,
    rabi_detuning_max_slope::Real,
)
    amp_times, amp_values = get_piecewise_linear_times_and_clocks(Ω, rabi_frequency_amplitude_max_slope)
    phase_times, phase_values = get_piecewise_linear_times_and_clocks(ϕ, rabi_frequency_phase_max_slope)
    detuning_times, detuning_values = get_piecewise_linear_times_and_clocks(Δ, rabi_detuning_max_slope)

    return EffectiveHamiltonian(;
        rydberg = RydbergHamiltonian(;
            rabi_frequency_amplitude = RydbergRabiFrequencyAmplitude(;
                global_value = RydbergRabiFrequencyAmplitudeGlobal(; times = amp_times, values = amp_values),
            ),
            rabi_frequency_phase = RydbergRabiFrequencyPhase(;
                global_value = RydbergRabiFrequencyPhaseGlobal(; times = phase_times, values = phase_values),
            ),
            detuning = RydbergDetuning(;
                global_value = RydbergDetuningGlobal(; times = detuning_times, values = detuning_values),
            ),
        ),
    )
end
