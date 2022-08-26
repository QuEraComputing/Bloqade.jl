


function get_rydberg_params(h::BloqadeExpr.RydbergHamiltonian)
    # extracts parameters from RydbergHamiltonian
    ϕ = nothing
    Ω = nothing
    Δ = nothing

    if h.rabi_term isa BloqadeExpr.SumOfX
        Ω = h.rabi_term.Ω.f
    elseif h.rabi_term isa BloqadeExpr.SumOfXPhase
        Ω = h.rabi_term.Ω.f
        ϕ = h.rabi_term.ϕ
    end

    if h.detuning_term isa BloqadeExpr.SumOfN
        Δ = h.detuning_term.Δ
    end

    return (h.rydberg_term.atoms,ϕ,Ω,Δ)
end

# parsing individual fields
# we split them up because each 
# field has different constraints. 

function parse_static_rydberg_Ω(param::Real,duration::Real,max_slope::Real,min_step::Real)
    step = max(min_step,abs(param/max_slope))

    clocks = Float64[0.0,step,duration-step,duration]
    values = Float64[0.0,param,param,0.0]
    waveform = piecewise_linear(;clocks,values)

    return (1.0,waveform)
end

function parse_static_rydberg_Ω(::Nothing,duration::Real,max_slope::Real,min_step::Real)

    clocks = Float64[0.0,duration]
    values = Float64[0.0,0.0]
    waveform = piecewise_linear(;clocks,values)

    return (1.0,waveform)
end

function parse_static_rydberg_Ω(::Vector{<:Real},duration::Real,max_slope::Real,min_step::Real)
    throw(ErrorException("Rabi field must be global drive."))
end

function parse_dynamic_rydberg_Ω(param::Waveform{F,T};duration=nothing) where {F,T<:Real}
    duration = (duration ≡ nothing ? param.duration : duration)
    
    if !(duration ≡ nothing) && param.duration != duration
        throw(ErrorException("Waveform durations do not match."))
    end

    if !isapprox(param(0.0), 0.0;atol=eps(),rtol=eps()) || !isapprox(param(duration),0.0;atol=eps(),rtol=eps())
        throw(ErrorException("Rabi Drive must start and end with value 0."))
    end

    return (1.0,param,duration)
end

function parse_dynamic_rydberg_Ω(::Vector{Waveform{F,T}};duration=nothing) where {F,T<:Real}
    throw(ErrorException("Rabi field must be global drive."))
end

parse_static_rydberg_Δ(param::Real,duration::Real,max_slope::Real,min_step::Real) = (1.0,piecewise_linear(;clocks=Float64[0.0,duration],values=Float64[param,param]))
parse_static_rydberg_Δ(::Nothing,duration::Real,max_slope::Real,min_step::Real) = (1.0,piecewise_linear(;clocks=Float64[0.0,duration],values=Float64[0.0,0.0]))
parse_static_rydberg_Δ(param::Vector{<:Real},duration::Real,max_slope::Real,min_step::Real) = (Float64[val for val in param],piecewise_linear(;clocks=Float64[0.0,duration],values=Float64[1,1]))

function parse_dynamic_rydberg_Δ(param::Waveform{F,T};duration=nothing) where {F,T<:Real}
    duration = (duration ≡ nothing ? param.duration : duration)
    
    if !(duration ≡ nothing) && param.duration != duration
        throw(ErrorException("Waveform durations do not match."))
    end

    return (1.0,param,duration)
end

function parse_dynamic_rydberg_Δ(param::Vector{Waveform{F,T}};duration=nothing) where {F,T<:Real}
    durations = [f.duration for f in param]
    
    duration = (duration ≡ nothing ? durations[1] : duration)

    if !(duration ≡ nothing) && !all(duration .== durations)
        throw(ErrorException("Waveform durations do not match."))
    end

    clock_samples = LinRange(0,duration,100)
    
    value_samples = zeros(length(param),length(clock_samples))

    for (i,f) in enumerate(param)
        for (j,clock) in enumerate(clock_samples)
            value_samples[i,j] = f(clock)
        end
    end
    # use SVD to determine if the waveforms are independent or not. 
    # if the there are more one nonzero singular value then there must be
    # multiple functions within the vector.
    u,s,vt = svd(value_samples)

    if any(s[2:end] .> s[1]*eps())
        throw(ErrorException("Local detuning waveforms cannot be decomposed into a product: Δ(i)⋅Δ(t)."))
    end 
    
    # use U vector to get the scal
    Amplitude = (u * Diagonal(s))[:,1]
    i = argmax(Amplitude)
    Amplitude ./= Amplitude[i]

    return (Amplitude,param[i],duration)
    
end
# ϕ has a combination of both Ω and Δ constraints. 
# e.g. has to be global, but can begin and end on non-zero values.

parse_static_rydberg_ϕ(param::Real,duration::Real,max_slope::Real,min_step::Real) = parse_static_rydberg_Δ(param,duration,max_slope,min_step)
parse_static_rydberg_ϕ(param::Nothing,duration::Real,max_slope::Real,min_step::Real) = parse_static_rydberg_Δ(param,duration,max_slope,min_step)
parse_static_rydberg_ϕ(param::Vector{<:Real},duration::Real,max_slope::Real,min_step::Real) = parse_static_rydberg_Ω(param,duration,max_slope,min_step)
parse_dynamic_rydberg_ϕ(param::Waveform{F,T};duration=nothing) where {F,T<:Real} = parse_dynamic_rydberg_Δ(param;duration)
parse_dynamic_rydberg_ϕ(param::Vector{Waveform{F,T}};duration=nothing) where {F,T<:Real} = parse_dynamic_rydberg_Ω(param;duration)


const ConstantParam = Union{Real,Nothing,Vector{<:Real}}
const DynamicParam = Union{Waveform{F,T} where {F,T<:Real},Vector{Waveform{F,T} where F} where T<:Real} 

# no dynamic parameeters must throw error because `duration` can't be determined
function parse_rydberg_fields(ϕ::ConstantParam,Ω::ConstantParam,Δ::ConstantParam) 
    throw(ErrorException("Schema requires at least one Waveform field."))
end
# one dynamic argument
function parse_rydberg_fields(ϕ::DynamicParam,Ω::ConstantParam,Δ::ConstantParam) end
function parse_rydberg_fields(ϕ::ConstantParam,Ω::DynamicParam,Δ::ConstantParam) end
function parse_rydberg_fields(ϕ::ConstantParam,Ω::ConstantParam,Δ::DynamicParam) end
# two dynamic arguments
function parse_rydberg_fields(ϕ::DynamicParam,Ω::DynamicParam,Δ::ConstantParam) end
function parse_rydberg_fields(ϕ::DynamicParam,Ω::ConstantParam,Δ::DynamicParam) end
function parse_rydberg_fields(ϕ::ConstantParam,Ω::DynamicParam,Δ::DynamicParam) end
# three dynamic arguments
function parse_rydberg_fields(ϕ::DynamicParam,Ω::DynamicParam,Δ::DynamicParam) end

parse_rydberg_fields(ϕ,Ω,Δ) = throw(UndefVarError("Unable to parse Rydberg coefficients for Schema conversion, please use Real/Nothing for constant coefficients and Waveforms dynamic coefficients."))
# function will extract parameters, check of waveforms are compatible with IR
# then descretize waveforms and return them to be parsed into 
# effective Hamiltonian. 

@inline convert_units(value::Real,from,to) = uconvert(to,Quantity(value,from))


function parse_analog_rydberg_params(h,params)
    atoms,ϕ,Ω,Δ = get_rydberg_params(h)
    # dispatch based on types
    ϕ,Ω,Δ,Δi = parse_rydberg_fields(ϕ,Ω,Δ)

    convert_time = x::Real -> convert_units(x,μs,s).val
    convert_rabi_amp = x::Real -> convert_units(x,rad*MHz,rad/s).val
    convert_rabi_phase = x::Real -> convert_units(x,rad,rad).val
    convert_detuning_local = x::Real -> convert_units(x,NoUnits,NoUnits).val

    # params are already unitful
    min_step = uconvert(μs,params.rabi_time_min_step).val
    

    max_slope = uconvert(rad*MHz/μs,params.rabi_frequency_phase_max_slope).val

    ϕ = BloqadeWaveforms.discretize(ϕ;
        max_slope=max_slope,
        min_step=min_step,
        tol=params.waveform_tolerance,
    )

    ϕ_clocks = [convert_time(clock) for clock in ϕ.f.clocks]
    ϕ_values = [convert_rabi_phase(value) for value in ϕ.f.values]
    ϕ_ir = piecewise_linear(;ϕ_clocks,ϕ_values)

    max_slope = uconvert(rad*MHz/μs,params.rabi_frequency_amplitude_max_slope).val
    
    Ω = BloqadeWaveforms.discretize(Ω;
        max_slope=max_slope,
        min_step=min_step,
        tol=params.waveform_tolerance,
    )

    Ω_clocks = [convert_time(clock) for clock in Ω.f.clocks]
    Ω_values = [convert_rabi_amp(value) for value in Ω.f.values]
    Ω_ir = piecewise_linear(;Ω_clocks,Ω_values)

    max_slope = uconvert(rad*MHz/μs,params.rabi_detuning_max_slope).val

    Δ = BloqadeWaveforms.discretize(Δ;
        max_slope=max_slope,
        min_step=min_step,
        tol=params.waveform_tolerance,
    )

    Δ_clocks = [convert_time(clock) for clock in Δ.f.clocks]
    Δ_values = [convert_rabi_amp(value) for value in Δ.f.values]
    Δ_ir = piecewise_linear(;Δ_clocks,Δ_values)
    Δi_ir = [convert_detuning_local(δ) for δ in Δi]

    return (atoms,Ω_ir,ϕ_ir,Δ_ir,Δi_ir)
end




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
    detuning_local = t.effective_hamiltonian.rydberg.detuning.local_value

    Ω = BloqadeWaveforms.piecewise_linear(; clocks=rabi_freq_amp.times, values=rabi_freq_amp.values)
    ϕ = BloqadeWaveforms.piecewise_linear(; clocks=rabi_freq_phase.times, values=rabi_freq_phase.values)
    Δ = BloqadeWaveforms.piecewise_linear(; clocks=detuning_global.times, values=detuning_global.values)

    if !isnothing(detuning_local)
        Δ_i = [δ_i*Δ for δ_i in detuning_local]
    else
        Δ_i = Δ
    end

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
to_json(h::BloqadeExpr.RydbergHamiltonian; kw...) = to_json(h, SchemaConversionParams(;kw...))
to_dict(h::BloqadeExpr.RydbergHamiltonian; kw...) = to_dict(h, SchemaConversionParams(;kw...))
to_schema(h::BloqadeExpr.RydbergHamiltonian; kw...) = to_schema(h, SchemaConversionParams(;kw...))

function to_json(h::BloqadeExpr.RydbergHamiltonian, params::SchemaConversionParams)
    return JSON.json(BloqadeSchema.to_dict(h, params))
end

function to_dict(h::BloqadeExpr.RydbergHamiltonian, params::SchemaConversionParams)
    return Configurations.to_dict(to_schema(h, params))
end

function to_schema(h::BloqadeExpr.RydbergHamiltonian, params::SchemaConversionParams)
    atoms,ϕ,Ω,Δ,Δ_i = parse_analog_rydberg_params(h,params)

    return TaskSpecification(;
        nshots=params.n_shots,
        lattice=to_lattice(atoms),
        effective_hamiltonian=to_hamiltonian(Ω, ϕ, Δ, Δ_i)
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
    Ω::Waveform{BloqadeWaveforms.PiecewiseLinear{T},T},
    ϕ::Waveform{BloqadeWaveforms.PiecewiseLinear{T},T},
    Δ::Waveform{BloqadeWaveforms.PiecewiseLinear{T},T},
    Δ_i::Vector{<:Real}) where {T<:Real}

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

function to_hamiltonian(
    Ω::Waveform{BloqadeWaveforms.PiecewiseLinear{T},T},
    ϕ::Waveform{BloqadeWaveforms.PiecewiseLinear{T},T},
    Δ::Waveform{BloqadeWaveforms.PiecewiseLinear{T},T},
    Δ_i::Real) where {T<:Real}

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
