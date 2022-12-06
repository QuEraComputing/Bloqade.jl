

@inline function convert_units(value::Real,from,to)

    val = uconvert(to,Quantity(value,from)).val
    return round(val;sigdigits=14) 
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

"""
    from_json(j::String)

Convert the JSON representation of a [`TaskSpecification`](@ref) instance to a 
[`TaskSpecification`](@ref)
"""
from_json(j::String) = BloqadeSchema.from_dict(JSON.parse(j))

"""
    from_dict(d::AbstractDict{String})

Convert the dictionary representation of a [`TaskSpecification`](@ref) instance, 
into a [`TaskSpecification`](@ref).
"""
function from_dict(d::AbstractDict{String})
    t = Configurations.from_dict(BloqadeSchema.TaskSpecification, d)
    return from_schema(t)
end

"""
    from_schema(t::TaskSpecification)

Converts `t` into valid `BloqadeExpr.RydbergHamiltonian` instance.
"""
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
        values=convert_units(rabi_freq_amp.values,rad/s,rad/μs)
    )
    ϕ = BloqadeWaveforms.piecewise_constant(;
        clocks=convert_units(rabi_freq_phase.times,s,μs),
        values=convert_units(rabi_freq_phase.values[1:end-1],rad,rad)
    )
    Δ = BloqadeWaveforms.piecewise_linear(;
        clocks=convert_units(detuning_global.times,s,μs),
        values=convert_units(detuning_global.values,rad/s,rad/μs)
    )
    if !isnothing(detuning_local)
        δ = BloqadeWaveforms.piecewise_linear(;
            clocks=convert_units(detuning_local.times,s,μs),
            values=convert_units(detuning_local.values,rad/s,rad/μs)
        )
        
        Δ_i = [Δ+δ_i*δ for (i,δ_i) in enumerate(detuning_local.lattice_site_coefficients) if t.lattice.filling[i] == 1]
    else
        Δ_i = Δ
    end

    return BloqadeExpr.rydberg_h(atoms; Δ=Δ_i, Ω=Ω, ϕ=ϕ)
end



"""
    to_json(h::AbstractBlock; kw...)
    to_json(h::BloqadeExpr.RydbergHamiltonian,params::SchemaTranslationParams)

Converts `h` and associated `params` into a JSON object.
If `params` is not explicitly provided as a `SchemaTranslationParams` instance, it is automatically built
from `nshots` and `device_capabilities`.

Validation is performed to ensure `h` is capable of being run on the machine. This can cause an
exception to be thrown should any violations be caught. Refer to Logs/Warnings/Exceptions below.

# Logs/Warnings/Exceptions

A `ValidationException` can be thrown which wraps a [`ValidationViolations`](@ref) instance.

`ValidationViolations` contains any constraint violations detected from [`to_schema`](@ref)

Violations include:

## Waveform Type
* ϕ is not of type `PiecewiseConstantWaveform`
* Ω and Δ are not of type `PiecewiseLinearWaveform`
## Atom Position
* Number of qubits requested exceeds what is supported by the device
* Atom positions exceed position resolution supported by the device
* The total width/height of the atom arrangement exceeds what is supported by the device
* The radial spacing between atoms is smaller than what is supported by the device
* The vertical row spacing between atoms is smaller than what is supported by the device
## General Waveform Constraints (apply to Ω, Δ, ϕ)
* duration exceeds device supported duration
* duration is smaller than device supported minimum time step
* smallest time step is smaller than supported smallest time step
* value is smaller than smallest supported value
* value is larger than largest supported value
## Ω Waveform specific constraints
* Slope exceeds largest supported slope
* Start and end values are not equal to 0.0 rad⋅MHz
## Δ Waveform specific constraints
* Slope exceeds largest supported slope
## ϕ Waveform specific constraints
* start value is not equal to 0.0 rad⋅MHz
## Miscellaneous Violations
* Number of shots is below minimum supported
* Number of shots exceeds maximum supported

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

"""
    to_dict(h::BloqadeExpr.RydbergHamiltonian; nshots::Int, device_capabilities::DeviceCapabilities=get_device_capabilities())
    to_dict(h::BloqadeExpr.RydbergHamiltonian,params::SchemaTranslationParams)

Converts `h` and associated `params` into the dictionary representation of a [`TaskSpecification`](@ref).
If `params` is not explicitly provided as a `SchemaTranslationParams` instance, it is automatically built
from `nshots` and `device_capabilities`.

Validation is performed to ensure `h` is capable of being run on the machine. This can cause an
exception to be thrown should any violations be caught. Refer to Logs/Warnings/Exceptions below.

# Logs/Warnings/Exceptions

A `ValidationException` can be thrown which wraps a [`ValidationViolations`](@ref) instance.

`ValidationViolations` contains any constraint violations detected from [`to_schema`](@ref)

Violations include:

## Waveform Type
* ϕ is not of type `PiecewiseConstantWaveform`
* Ω and Δ are not of type `PiecewiseLinearWaveform`
## Atom Position
* Number of qubits requested exceeds what is supported by the device
* Atom positions exceed position resolution supported by the device
* The total width/height of the atom arrangement exceeds what is supported by the device
* The radial spacing between atoms is smaller than what is supported by the device
* The vertical row spacing between atoms is smaller than what is supported by the device
## General Waveform Constraints (apply to Ω, Δ, ϕ)
* duration exceeds device supported duration
* duration is smaller than device supported minimum time step
* smallest time step is smaller than supported smallest time step
* value is smaller than smallest supported value
* value is larger than largest supported value
## Ω Waveform specific constraints
* Slope exceeds largest supported slope
* Start and end values are not equal to 0.0 rad⋅MHz
## Δ Waveform specific constraints
* Slope exceeds largest supported slope
## ϕ Waveform specific constraints
* start value is not equal to 0.0 rad⋅MHz
## Miscellaneous Violations
* Number of shots is below minimum supported
* Number of shots exceeds maximum supported

"""
to_dict(h::BloqadeExpr.RydbergHamiltonian; kw...) = to_dict(h,SchemaTranslationParams(;kw...))

"""
    to_schema(h::BloqadeExpr.RydbergHamiltonian; nshots::Int, device_capabilities::DeviceCapabilities=get_device_capabilities())
    to_schema(h::BloqadeExpr.RydbergHamiltonian, params::SchemaTranslationParams)

Converts `h` to a `TaskSpecification` instance with `params`. If params is not explicitly constructed, 
it will be built automatically from `nshots` and `device_capabilities`. 

Validation is performed to ensure `h` is capable of being run on the machine. This can cause an
exception to be thrown should any violations be caught. Refer to Logs/Warnings/Exceptions below.

# Logs/Warnings/Exceptions

If any violations of `device_capabilities` are detected, a `ValidationException` is thrown which wraps
a [`ValidationViolations`](@ref) instance.

Violations include:

## Waveform Type
* ϕ is not of type `PiecewiseConstantWaveform`
* Ω and Δ are not of type `PiecewiseLinearWaveform`
## Atom Position
* Number of qubits requested exceeds what is supported by the device
* Atom positions exceed position resolution supported by the device
* The total width/height of the atom arrangement exceeds what is supported by the device
* The vertical row spacing between atoms is smaller than what is supported by the device
* The radial spacing between atoms is smaller than what is supported by the device
## General Waveform Constraints (apply to Ω, Δ, ϕ)
* duration exceeds device supported duration
* duration is smaller than device supported minimum time step
* smallest time step is smaller than supported smallest time step
* value is smaller than smallest supported value
* value is larger than largest supported value
## Ω Waveform specific constraints
* Slope exceeds largest supported slope
* Start and end values are not equal to 0.0 rad⋅MHz
## Δ Waveform specific constraints
* Slope exceeds largest supported slope
## ϕ Waveform specific constraints
* start value is not equal to 0.0 rad⋅MHz
## Miscellaneous Violations
* Number of shots is below minimum supported
* Number of shots exceeds maximum supported
"""
to_schema(h::BloqadeExpr.RydbergHamiltonian; kw...) = to_schema(h,SchemaTranslationParams(;kw...))

function to_json(h::BloqadeExpr.RydbergHamiltonian,params::SchemaTranslationParams)
    return JSON.json(BloqadeSchema.to_dict(h,params))
end

function to_dict(h::BloqadeExpr.RydbergHamiltonian,params::SchemaTranslationParams)
    return Configurations.to_dict(to_schema(h,params))
end

function to_schema(h::BloqadeExpr.RydbergHamiltonian, params::SchemaTranslationParams)
    atoms,ϕ,Ω,Δ,δ,Δi = schema_parse_rydberg_fields(h)

    violations = validate_analog_params(atoms,ϕ,Ω,Δ,δ,Δi,params.device_capabilities)

    params.n_shots < params.device_capabilities.task.number_shots_min && push!(violations.misc_violations,
        "n_shots $(params.n_shots) is less than minimum value $(params.device_capabilities.task.number_shots_min)"
    )
    params.n_shots > params.device_capabilities.task.number_shots_max && push!(violations.misc_violations,
        "n_shots $(params.n_shots) is exceeds maximum value $(params.device_capabilities.task.number_shots_max)"
    )

    if !isempty(violations)
        throw(ValidationException(violations))
    end

    ϕ =(
        clocks=ϕ.f.clocks,
        values=[ϕ.f.values...,ϕ.f.values[end]] # add extra element at the end 
    )

    Ω = (
        clocks=Ω.f.clocks,
        values=Ω.f.values
    )

    Δ = (
        clocks=Δ.f.clocks,
        values=Δ.f.values
    )

    δ = if !isnothing(δ) 
        (
            clocks=δ.f.clocks,
            values=δ.f.values
        )
    end

    return TaskSpecification(;
        nshots=params.n_shots,
        lattice=to_lattice(atoms),
        effective_hamiltonian=to_hamiltonian(Ω, ϕ, Δ, δ, Δi)
    )
end

"""
    to_json_no_validation(lattice::Union{Vector,Lattice};
        ϕ::Maybe{PiecewiseConstantWaveform}=nothing,
        Ω::Maybe{PiecewiseLinearWaveform}=nothing,
        Δ::Maybe{PiecewiseLinearWaveform}=nothing,
        δ::Maybe{PiecewiseLinearWaveform}=nothing,
        Δi::Maybe{Vector{Number}}=nothing,kw...)

Converts `lattice`, `ϕ`, `Δ`, `δ`, and `Δi` to a JSON representation of a `TaskSpecification` instance
WITHOUT ensuring the provided values are capable of being executed on the machine (fit within the 
constraints of the device's capabilities)

See also [`to_json`](@ref)
"""
function to_json_no_validation(lattice::Union{Vector,Lattice};
    ϕ::Maybe{PiecewiseConstantWaveform}=nothing,
    Ω::Maybe{PiecewiseLinearWaveform}=nothing,
    Δ::Maybe{PiecewiseLinearWaveform}=nothing,
    δ::Maybe{PiecewiseLinearWaveform}=nothing,
    Δi::Maybe{Vector{Number}}=nothing,kw...)
    schema = to_schema_no_validation(lattice,ϕ,Ω,Δ,δ,Δi,SchemaTranslationParams(;kw...))
    return JSON.json(Configurations.to_dict(schema))
end

"""
    to_schema_no_validation(lattice::Union{Vector,Lattice};
        ϕ::Maybe{PiecewiseConstantWaveform}=nothing,
        Ω::Maybe{PiecewiseLinearWaveform}=nothing,
        Δ::Maybe{PiecewiseLinearWaveform}=nothing,
        δ::Maybe{PiecewiseLinearWaveform}=nothing,
        Δi::Maybe{Vector{Number}}=nothing, 
        nshots::Int,
        device_capabilities::DeviceCapabilities=get_device_capabilities())
    to_schema_no_validation(lattice::Union{Vector,Lattice},
        ϕ::Maybe{PiecewiseConstantWaveform},
        Ω::Maybe{PiecewiseLinearWaveform},
        Δ::Maybe{PiecewiseLinearWaveform},
        δ::Maybe{PiecewiseLinearWaveform},
        Δi::Maybe{Vector{Number}}, 
        params::SchemaTranslationParams)

Converts `lattice`, `ϕ`, `Δ`, `δ`, and `Δi` to a `TaskSpecification` instance
WITHOUT ensuring the provided values are capable of being executed on the machine (fit within the 
constraints of the device's capabilities).

If `params` is not already provided, it is constructed automatically from `nshots::Int`
and `device_capabilities`.

See also [`to_schema`](@ref)
"""
function to_schema_no_validation(lattice::Union{Vector,Lattice};
    ϕ::Maybe{PiecewiseConstantWaveform}=nothing,
    Ω::Maybe{PiecewiseLinearWaveform}=nothing,
    Δ::Maybe{PiecewiseLinearWaveform}=nothing,
    δ::Maybe{PiecewiseLinearWaveform}=nothing,
    Δi::Maybe{Vector{Number}}=nothing, 
    kw...)
    return to_schema_no_validation(lattice,ϕ,Ω,Δ,δ,Δi,SchemaTranslationParams(;kw...))
end

function to_schema_no_validation(lattice::Union{Vector,Lattice},
    ϕ::Maybe{PiecewiseConstantWaveform},
    Ω::Maybe{PiecewiseLinearWaveform},
    Δ::Maybe{PiecewiseLinearWaveform},
    δ::Maybe{PiecewiseLinearWaveform},
    Δi::Maybe{Vector{Number}}, 
    params::SchemaTranslationParams)
    
    duration = params.device_capabilities.rydberg.global_value.time_delta_min
    # take maximum value for duration
    
    for wf in [Ω,Δ,ϕ]
        duration =  (!isnothing(wf) ? max(wf.duration,duration) : duration)
    end

    ϕ = if !isnothing(ϕ)
        (
            clocks=ϕ.f.clocks,
            values=[ϕ.f.values...,ϕ.f.values[end]] # add extra element at the end 
        )
    else
        (
            clocks = [0.0,0.0],
            values = [0.0,duration]
        )
    end

    Ω = if !isnothing(Ω)
        (
            clocks=Ω.f.clocks,
            values=Ω.f.values
        )
    else
        (
            clocks = [0.0,0.0],
            values = [0.0,duration]
        )
    end

    Δ =  if !isnothing(Δ) 
        (
            clocks=Δ.f.clocks,
            values=Δ.f.values
        )
    else
        (
            clocks = [0.0,0.0],
            values = [0.0,duration]
        )
    end

    δ = if !isnothing(δ)  
        (
            clocks=δ.f.clocks,
            values=δ.f.values
        )
    end

    return TaskSpecification(;
        nshots=params.n_shots,
        lattice=to_lattice(lattice),
        effective_hamiltonian=to_hamiltonian(Ω, ϕ, Δ, δ, Δi)
    )

end

to_lattice(lattice::Lattice) = Lattice(
    sites = [convert_units.(site,μm,m) for site in lattice.sites],
    filling = lattice.filling
)

function to_lattice(atoms::Vector)
    coords = map(atoms) do coord
        length(coord) == 1 && return convert_units.((coord[1], 0),μm,m)
        return convert_units.(coord,μm,m)
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
                global_value = RydbergRabiFrequencyAmplitudeGlobal(; 
                    times = convert_units.(Ω.clocks,μs,s), 
                    values = convert_units.(Ω.values,rad*MHz,rad/s)
                ),
            ),
            rabi_frequency_phase = RydbergRabiFrequencyPhase(;
                global_value = RydbergRabiFrequencyPhaseGlobal(; 
                    times = convert_units.(ϕ.clocks,μs,s),
                    values = convert_units.(ϕ.values,rad,rad)
                ),
            ),
            detuning = RydbergDetuning(;
                global_value = RydbergDetuningGlobal(; 
                    times = convert_units.(Δ.clocks,μs,s), 
                    values = convert_units.(Δ.values,rad*MHz,rad/s)
                ),
                local_value = RydbergDetuningLocal(; 
                    times = convert_units.(δ.clocks,μs,s), 
                    values = convert_units.(δ.values,rad*MHz,rad/s), 
                    lattice_site_coefficients=Δ_i
                )
            ),
        ),
    )
end

function to_hamiltonian(
    Ω::NamedTuple,
    ϕ::NamedTuple,
    Δ::NamedTuple,
    δ::Nothing,
    Δ_i::Maybe{Real})

    return EffectiveHamiltonian(;
        rydberg = RydbergHamiltonian(;
            rabi_frequency_amplitude = RydbergRabiFrequencyAmplitude(;
                global_value = RydbergRabiFrequencyAmplitudeGlobal(; 
                    times = convert_units.(Ω.clocks,μs,s), 
                    values = convert_units.(Ω.values,rad/μs,rad/s)
                ),
            ),
            rabi_frequency_phase = RydbergRabiFrequencyPhase(;
                global_value = RydbergRabiFrequencyPhaseGlobal(; 
                    times = convert_units.(ϕ.clocks,μs,s),
                    values = convert_units.(ϕ.values,rad,rad)
                ),
            ),
            detuning = RydbergDetuning(;
                global_value = RydbergDetuningGlobal(; 
                    times = convert_units.(Δ.clocks,μs,s), 
                    values = convert_units.(Δ.values,rad/μs,rad/s)
                )
            ),
        ),
    )
end
