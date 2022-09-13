

@inline function convert_units(value::Real,from,to)

    val = uconvert(to,Quantity(value,from)).val
    return round(val;sigdigits=15) 
end

function convert_units(x::AbstractArray{S},from,to) where {S<:Real}
    y = similar(x)
    @inbounds for i in eachindex(x)
        y[i] = convert_units(x[i],from,to)
    end
    return y
end



"""
    execute(j::String)

Executes a task given as a JSON string in the task specification API format, and returns a JSON string of the result
"""

function execute(j::String)
    execute(JSON.parse(j))
end

"""
    execute(j::Dict)

Executes a task given as a Dict in the task specification API format, and returns a JSON string of the result
"""
function execute(dict::AbstractDict{String})
    execute(Configurations.from_dict(BloqadeSchema.TaskSpecification,dict))
end


"""
    execute(j::TaskSpecification)

Executes a task given as a TaskSpecification object in the task specification API format, and returns a JSON string of the result
"""
function execute(task::TaskSpecification)
    h = from_schema(task)
        
    atoms,ϕ,Ω,Δ = get_rydberg_params(h)
    total_time = ϕ.duration
    n_atoms = length(atoms)

    reg = zero_state(n_atoms)
    problem = SchrodingerProblem(reg, total_time, h)
    BloqadeExpr.emulate!(problem)
    bitstrings = reg |> measure(; nshots=task.nshots)

    return JSON.json(Configurations.to_dict(to_task_output(bitstrings)))
end




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
    atoms = (site for (i,site) in enumerate(t.lattice.sites) if t.lattice.filling[i] == 1)

    atoms = map(atoms) do pos 
        return convert_units.(pos,m,μm)
    end

    rabi_freq_amp = t.effective_hamiltonian.rydberg.rabi_frequency_amplitude.global_value
    rabi_freq_phase = t.effective_hamiltonian.rydberg.rabi_frequency_phase.global_value
    detuning_global = t.effective_hamiltonian.rydberg.detuning.global_value
    detuning_local = t.effective_hamiltonian.rydberg.detuning.local_value

    Ω = BloqadeWaveforms.piecewise_linear(; 
        clocks=convert_units(rabi_freq_amp.times,s,μs), 
        values=convert_units(rabi_freq_amp.values,rad/s,rad*MHz)
    )
    ϕ = BloqadeWaveforms.piecewise_linear(;
        clocks=convert_units(rabi_freq_phase.times,s,μs),
        values=convert_units(rabi_freq_phase.values,rad,rad)
    )
    Δ = BloqadeWaveforms.piecewise_linear(;
        clocks=convert_units(detuning_global.times,s,μs),
        values=convert_units(detuning_global.values,rad/s,rad*MHz)
    )
    if !isnothing(detuning_local)
        δ = BloqadeWaveforms.piecewise_linear(;
            clocks=convert_units(detuning_local.times,s,μs),
            values=convert_units(detuning_local.values,rad/s,rad*MHz)
        )
        
        Δ_i = [Δ+δ_i*δ for (i,δ_i) in enumerate(detuning_local.lattice_site_coefficients) if t.lattice.filling[i] == 1]
    else
        Δ_i = Δ
    end

    all(norm.(Δ) .< eps()) && Δ = nothing
    norm(Ω) < eps() && Ω = nothing
    norm(ϕ) < eps() && ϕ = nothing

    return BloqadeExpr.rydberg_h(atoms; Δ=Δ_i, Ω=Ω, ϕ=ϕ)
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
to_json(h::BloqadeExpr.RydbergHamiltonian; kw...) = to_json(h,SchemaTranslationParams(;kw...))
to_dict(h::BloqadeExpr.RydbergHamiltonian; kw...) = to_dict(h,SchemaTranslationParams(;kw...))
to_schema(h::BloqadeExpr.RydbergHamiltonian; kw...) = to_schema(h,SchemaTranslationParams(;kw...))

function to_json(h::BloqadeExpr.RydbergHamiltonian,params::SchemaTranslationParams)
    return JSON.json(BloqadeSchema.to_dict(h,params))
end

function to_dict(h::BloqadeExpr.RydbergHamiltonian,params::SchemaTranslationParams)
    return Configurations.to_dict(to_schema(h,params))
end



function to_schema(h::BloqadeExpr.RydbergHamiltonian, params::SchemaTranslationParams)
    atoms,ϕ,Ω,Δ,info = hardware_transform_parse(h;params.device_capabilities)

    # extract Detuning mask
    Δ = info.Δ_mask.Δ
    δ = info.Δ_mask.δ
    Δi = info.Δ_mask.Δi

    if params.transform_info
        ϕ_err = info.ϕ
        Ω_err = info.Ω
        Δ_err = info.Δ
        @info "Hardware transform report: ∫dt |ϕ(t)-ϕ_hw(t)| = $ϕ_err"
        @info "Hardware transform report: ∫dt |Ω(t)-Ω_hw(t)| = $Ω_err"
        @info "Hardware transform report: ∫dt |Δ(t)-Δ_hw(t)| = $Δ_err"
    end

    # validate components to save conversions
    validate_lattice(atoms,params.warn,params.device_capabilities)
    rydberg_capabilities = get_rydberg_capabilities(;device_capabilities=params.device_capabilities)

    validate_ϕ(ϕ,params.warn,rydberg_capabilities.ϕ)
    validate_Ω(Ω,params.warn,rydberg_capabilities.Ω)
    validate_Δ(Δ,params.warn,rydberg_capabilities.Δ)
    if !isnothing(δ)
        validate_δ(δ,Δi,warn,rydberg_capabilities.δ)
    end        
    atoms = map(atoms) do pos 
        return convert_units.(pos,μm,m)
    end
    
    ϕ_clocks = convert_units(ϕ.f.clocks,μs,s) 
    ϕ_values = convert_units(ϕ.f.values,rad,rad)
    Ω_clocks = convert_units(Ω.f.clocks,μs,s)
    Ω_values = convert_units(Ω.f.values,rad*MHz,rad/s)
    Δ_clocks = convert_units(Δ.f.clocks,μs,s)
    Δ_values = convert_units(Δ.f.values,rad*MHz,rad/s)

    ϕ =(
        clocks=ϕ_clocks,
        values=ϕ_values
    )


    Ω = (
        clocks=Ω_clocks,
        values=Ω_values
    )


    Δ = (
        clocks=Δ_clocks,
        values=Δ_values
    )

    if !isnothing(δ)
        δ_clocks = convert_units(δ.f.clocks,μs,s)
        δ_values = convert_units(δ.f.values,rad*MHz,rad/s)
        
        δ = (
            clocks=δ_clocks,
            values=δ_values
        )
    else

    end

    return TaskSpecification(;
        nshots=params.n_shots,
        lattice=to_lattice(atoms),
        effective_hamiltonian=to_hamiltonian(Ω, ϕ, Δ, δ, Δi)
    )
end


function to_lattice(atoms::Vector)
    coords = map(atoms) do coord
        length(coord) == 1 && return (coord[1], 0)
        return coord
    end
    return Lattice(; sites = coords, filling = vec(ones(length(coords), 1)))
end



function to_hamiltonian(
    Ω::NamedTuple,
    ϕ::NamedTuple,
    Δ::NamedTuple,
    δ::NamedTuple,
    Δ_i::Vector{<:Real}) 

    return EffectiveHamiltonian(;
        rydberg = RydbergHamiltonian(;
            rabi_frequency_amplitude = RydbergRabiFrequencyAmplitude(;
                global_value = RydbergRabiFrequencyAmplitudeGlobal(; times = Ω.clocks, values = Ω.values),
            ),
            rabi_frequency_phase = RydbergRabiFrequencyPhase(;
                global_value = RydbergRabiFrequencyPhaseGlobal(; times = ϕ.clocks, values = ϕ.values),
            ),
            detuning = RydbergDetuning(;
                global_value = RydbergDetuningGlobal(; times = Δ.clocks, values = Δ.values),
                local_value = RydbergDetuningLocal(; times = δ.clocks, values = δ.values, lattice_site_coefficients=Δ_i)
            ),
        ),
    )
end

function to_hamiltonian(
    Ω::NamedTuple,
    ϕ::NamedTuple,
    Δ::NamedTuple,
    δ::Nothing,
    Δ_i::Real)

    return EffectiveHamiltonian(;
        rydberg = RydbergHamiltonian(;
            rabi_frequency_amplitude = RydbergRabiFrequencyAmplitude(;
                global_value = RydbergRabiFrequencyAmplitudeGlobal(; times = Ω.clocks, values = Ω.values),
            ),
            rabi_frequency_phase = RydbergRabiFrequencyPhase(;
                global_value = RydbergRabiFrequencyPhaseGlobal(; times = ϕ.clocks, values = ϕ.values),
            ),
            detuning = RydbergDetuning(;
                global_value = RydbergDetuningGlobal(; times = Δ.clocks, values = Δ.values),
                # local_value = RydbergDetuningLocal(; times = δ.clocks, values = δ.values, lattice_site_coefficients=Δ_i)
            ),
        ),
    )
end
