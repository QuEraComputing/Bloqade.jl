

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

    return waveform
end

function parse_static_rydberg_Ω(::Nothing,duration::Real,max_slope::Real,min_step::Real)

    clocks = Float64[0.0,duration]
    values = Float64[0.0,0.0]
    waveform = piecewise_linear(;clocks,values)

    return waveform
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

@inline convert_units(value::Real,from,to) = uconvert(to,Quantity(value,from))

# no dynamic parameeters must throw error because `duration` can't be determined
function parse_analog_rydberg_fields(ϕ::ConstantParam,Ω::ConstantParam,Δ::ConstantParam,params) 
    throw(ErrorException("Schema requires at least one Waveform field."))
end

# one dynamic argument
function parse_analog_rydberg_fields(ϕ::DynamicParam,Ω::ConstantParam,Δ::ConstantParam,params) 
    convert_time = x::Real -> convert_units(x,μs,s).val
    convert_rabi_amp = x::Real -> convert_units(x,rad*MHz,rad/s).val
    convert_rabi_phase = x::Real -> convert_units(x,rad,rad).val

    # params are already unitful
    min_step = uconvert(μs,params.rabi_time_min_step).val
    ϕ_max_slope = uconvert(rad*MHz/μs,params.rabi_frequency_phase_max_slope).val
    Ω_max_slope = uconvert(rad*MHz/μs,params.rabi_frequency_amplitude_max_slope).val
    Δ_max_slope = uconvert(rad*MHz/μs,params.rabi_detuning_max_slope).val

    ϕ,duration = parse_dynamic_rydberg_ϕ(ϕ)
    Ω = parse_static_rydberg_Ω(Ω,duration,Ω_max_slope,min_step)
    Δ_local,Δ = parse_static_rydberg_Δ(Δ,duration,Δ_max_slope,min_step)

    ϕ = BloqadeWaveforms.discretize(ϕ;
        max_slope=ϕ_max_slope,
        min_step=min_step,
        tol=params.waveform_tolerance,
    )

    ϕ = piecewise_linear(;
        clocks=[convert_time(clock) for clock in ϕ.f.clocks],
        values=[convert_rabi_phase(value) for value in ϕ.f.values]
    )

    Ω = piecewise_linear(;
        clocks=[convert_time(clock) for clock in Ω.f.clocks],
        values=[convert_rabi_amp(value) for value in Ω.f.values]
    )

    Δ = piecewise_linear(;
        clocks=[convert_time(clock) for clock in Δ.f.clocks],
        values=[convert_rabi_amp(value) for value in Δ.f.values]
    )
    return (ϕ,Ω,Δ,Δ_local)
end

function parse_analog_rydberg_fields(ϕ::ConstantParam,Ω::DynamicParam,Δ::ConstantParam,params)
    convert_time = x::Real -> convert_units(x,μs,s).val
    convert_rabi_amp = x::Real -> convert_units(x,rad*MHz,rad/s).val
    convert_rabi_phase = x::Real -> convert_units(x,rad,rad).val

    # params are already unitful
    min_step = uconvert(μs,params.rabi_time_min_step).val
    ϕ_max_slope = uconvert(rad*MHz/μs,params.rabi_frequency_phase_max_slope).val
    Ω_max_slope = uconvert(rad*MHz/μs,params.rabi_frequency_amplitude_max_slope).val
    Δ_max_slope = uconvert(rad*MHz/μs,params.rabi_detuning_max_slope).val

    Ω,duration = parse_dynamic_rydberg_ϕ(Ω)
    ϕ = parse_static_rydberg_ϕ(ϕ,duration,ϕ_max_slope,min_step)
    Δ_local,Δ = parse_static_rydberg_Δ(Δ,duration,Δ_max_slope,min_step)

    Ω = BloqadeWaveforms.discretize(Ω;
        max_slope=Ω_max_slope,
        min_step=min_step,
        tol=params.waveform_tolerance,
    )

    ϕ = piecewise_linear(;
        clocks=[convert_time(clock) for clock in ϕ.f.clocks],
        values=[convert_rabi_phase(value) for value in ϕ.f.values]
    )

    Ω = piecewise_linear(;
        clocks=[convert_time(clock) for clock in Ω.f.clocks],
        values=[convert_rabi_amp(value) for value in Ω.f.values]
    )

    Δ = piecewise_linear(;
        clocks=[convert_time(clock) for clock in Δ.f.clocks],
        values=[convert_rabi_amp(value) for value in Δ.f.values]
    )
    return (ϕ,Ω,Δ,Δ_local)
end

function parse_analog_rydberg_fields(ϕ::ConstantParam,Ω::ConstantParam,Δ::DynamicParam,params) 
    convert_time = x::Real -> convert_units(x,μs,s).val
    convert_rabi_amp = x::Real -> convert_units(x,rad*MHz,rad/s).val
    convert_rabi_phase = x::Real -> convert_units(x,rad,rad).val

    # params are already unitful
    min_step = uconvert(μs,params.rabi_time_min_step).val
    ϕ_max_slope = uconvert(rad*MHz/μs,params.rabi_frequency_phase_max_slope).val
    Ω_max_slope = uconvert(rad*MHz/μs,params.rabi_frequency_amplitude_max_slope).val
    Δ_max_slope = uconvert(rad*MHz/μs,params.rabi_detuning_max_slope).val

    Δ_local,Δ,duration = parse_dynamic_rydberg_Δ(Δ)
    ϕ = parse_static_rydberg_ϕ(ϕ,duration,ϕ_max_slope,min_step)
    Ω = parse_static_rydberg_Ω(Ω,duration,Ω_max_slope,min_step)

    Δ = BloqadeWaveforms.discretize(Δ;
        max_slope=Δ_max_slope,
        min_step=min_step,
        tol=params.waveform_tolerance,
    )

    ϕ = piecewise_linear(;
        clocks=[convert_time(clock) for clock in ϕ.f.clocks],
        values=[convert_rabi_phase(value) for value in ϕ.f.values]
    )

    Ω = piecewise_linear(;
        clocks=[convert_time(clock) for clock in Ω.f.clocks],
        values=[convert_rabi_amp(value) for value in Ω.f.values]
    )

    Δ = piecewise_linear(;
        clocks=[convert_time(clock) for clock in Δ.f.clocks],
        values=[convert_rabi_amp(value) for value in Δ.f.values]
    )
    return (ϕ,Ω,Δ,Δ_local)
end

# two dynamic arguments
function parse_analog_rydberg_fields(ϕ::DynamicParam,Ω::DynamicParam,Δ::ConstantParam,params) 
    convert_time = x::Real -> convert_units(x,μs,s).val
    convert_rabi_amp = x::Real -> convert_units(x,rad*MHz,rad/s).val
    convert_rabi_phase = x::Real -> convert_units(x,rad,rad).val

    # params are already unitful
    min_step = uconvert(μs,params.rabi_time_min_step).val
    ϕ_max_slope = uconvert(rad*MHz/μs,params.rabi_frequency_phase_max_slope).val
    Ω_max_slope = uconvert(rad*MHz/μs,params.rabi_frequency_amplitude_max_slope).val
    Δ_max_slope = uconvert(rad*MHz/μs,params.rabi_detuning_max_slope).val

    ϕ,duration = parse_dynamic_rydberg_ϕ(ϕ)
    Ω,duration = parse_dynamic_rydberg_Ω(Ω;duration)
    Δ_local,Δ = parse_static_rydberg_Δ(Δ,duration,Δ_max_slope,min_step)

    ϕ = BloqadeWaveforms.discretize(ϕ;
        max_slope=ϕ_max_slope,
        min_step=min_step,
        tol=params.waveform_tolerance,
    )

    Ω = BloqadeWaveforms.discretize(Ω;
        max_slope=Ω_max_slope,
        min_step=min_step,
        tol=params.waveform_tolerance,
    )

    ϕ = piecewise_linear(;
        clocks=[convert_time(clock) for clock in ϕ.f.clocks],
        values=[convert_rabi_phase(value) for value in ϕ.f.values]
    )

    Ω = piecewise_linear(;
        clocks=[convert_time(clock) for clock in Ω.f.clocks],
        values=[convert_rabi_amp(value) for value in Ω.f.values]
    )

    Δ = piecewise_linear(;
        clocks=[convert_time(clock) for clock in Δ.f.clocks],
        values=[convert_rabi_amp(value) for value in Δ.f.values]
    )
    return (ϕ,Ω,Δ,Δ_local)
end

function parse_analog_rydberg_fields(ϕ::DynamicParam,Ω::ConstantParam,Δ::DynamicParam,params) 
    convert_time = x::Real -> convert_units(x,μs,s).val
    convert_rabi_amp = x::Real -> convert_units(x,rad*MHz,rad/s).val
    convert_rabi_phase = x::Real -> convert_units(x,rad,rad).val

    # params are already unitful
    min_step = uconvert(μs,params.rabi_time_min_step).val
    ϕ_max_slope = uconvert(rad*MHz/μs,params.rabi_frequency_phase_max_slope).val
    Ω_max_slope = uconvert(rad*MHz/μs,params.rabi_frequency_amplitude_max_slope).val
    Δ_max_slope = uconvert(rad*MHz/μs,params.rabi_detuning_max_slope).val

    ϕ,duration = parse_dynamic_rydberg_ϕ(ϕ)
    Δ_local,Δ,duration = parse_dynamic_rydberg_Δ(Δ;duration)
    Ω = parse_static_rydberg_Ω(Ω,duration,Ω_max_slope,min_step)

    ϕ = BloqadeWaveforms.discretize(ϕ;
        max_slope=ϕ_max_slope,
        min_step=min_step,
        tol=params.waveform_tolerance,
    )

    Δ = BloqadeWaveforms.discretize(Δ;
        max_slope=Δ_max_slope,
        min_step=min_step,
        tol=params.waveform_tolerance,
    )

    ϕ = piecewise_linear(;
        clocks=[convert_time(clock) for clock in ϕ.f.clocks],
        values=[convert_rabi_phase(value) for value in ϕ.f.values]
    )

    Ω = piecewise_linear(;
        clocks=[convert_time(clock) for clock in Ω.f.clocks],
        values=[convert_rabi_amp(value) for value in Ω.f.values]
    )

    Δ = piecewise_linear(;
        clocks=[convert_time(clock) for clock in Δ.f.clocks],
        values=[convert_rabi_amp(value) for value in Δ.f.values]
    )
    return (ϕ,Ω,Δ,Δ_local)
end

function parse_analog_rydberg_fields(ϕ::ConstantParam,Ω::DynamicParam,Δ::DynamicParam,params) 
    convert_time = x::Real -> convert_units(x,μs,s).val
    convert_rabi_amp = x::Real -> convert_units(x,rad*MHz,rad/s).val
    convert_rabi_phase = x::Real -> convert_units(x,rad,rad).val

    # params are already unitful
    min_step = uconvert(μs,params.rabi_time_min_step).val
    ϕ_max_slope = uconvert(rad*MHz/μs,params.rabi_frequency_phase_max_slope).val
    Ω_max_slope = uconvert(rad*MHz/μs,params.rabi_frequency_amplitude_max_slope).val
    Δ_max_slope = uconvert(rad*MHz/μs,params.rabi_detuning_max_slope).val

    Ω,duration = parse_dynamic_rydberg_Ω(Ω)
    Δ_local,Δ,duration = parse_dynamic_rydberg_Δ(Δ;duration)
    ϕ = parse_static_rydberg_ϕ(ϕ,duration,ϕ_max_slope,min_step)

    Ω = BloqadeWaveforms.discretize(Ω;
        max_slope=Ω_max_slope,
        min_step=min_step,
        tol=params.waveform_tolerance,
    )

    Δ = BloqadeWaveforms.discretize(Δ;
        max_slope=Δ_max_slope,
        min_step=min_step,
        tol=params.waveform_tolerance,
    )

    ϕ = piecewise_linear(;
        clocks=[convert_time(clock) for clock in ϕ.f.clocks],
        values=[convert_rabi_phase(value) for value in ϕ.f.values]
    )

    Ω = piecewise_linear(;
        clocks=[convert_time(clock) for clock in Ω.f.clocks],
        values=[convert_rabi_amp(value) for value in Ω.f.values]
    )

    Δ = piecewise_linear(;
        clocks=[convert_time(clock) for clock in Δ.f.clocks],
        values=[convert_rabi_amp(value) for value in Δ.f.values]
    )
    return (ϕ,Ω,Δ,Δ_local)
end

# three dynamic arguments
function parse_analog_rydberg_fields(ϕ::DynamicParam,Ω::DynamicParam,Δ::DynamicParam,params) 
    convert_time = x::Real -> convert_units(x,μs,s).val
    convert_rabi_amp = x::Real -> convert_units(x,rad*MHz,rad/s).val
    convert_rabi_phase = x::Real -> convert_units(x,rad,rad).val

    # params are already unitful
    min_step = uconvert(μs,params.rabi_time_min_step).val
    ϕ_max_slope = uconvert(rad*MHz/μs,params.rabi_frequency_phase_max_slope).val
    Ω_max_slope = uconvert(rad*MHz/μs,params.rabi_frequency_amplitude_max_slope).val
    Δ_max_slope = uconvert(rad*MHz/μs,params.rabi_detuning_max_slope).val

    ϕ,duration = parse_dynamic_rydberg_Ω(ϕ)
    Ω,duration = parse_dynamic_rydberg_Ω(Ω;duration)
    Δ_local,Δ,duration = parse_dynamic_rydberg_Δ(Δ;duration)
    
    ϕ = BloqadeWaveforms.discretize(ϕ;
        max_slope=ϕ_max_slope,
        min_step=min_step,
        tol=params.waveform_tolerance,
    )

    Ω = BloqadeWaveforms.discretize(Ω;
        max_slope=Ω_max_slope,
        min_step=min_step,
        tol=params.waveform_tolerance,
    )

    Δ = BloqadeWaveforms.discretize(Δ;
        max_slope=Δ_max_slope,
        min_step=min_step,
        tol=params.waveform_tolerance,
    )

    ϕ = piecewise_linear(;
        clocks=[convert_time(clock) for clock in ϕ.f.clocks],
        values=[convert_rabi_phase(value) for value in ϕ.f.values]
    )

    Ω = piecewise_linear(;
        clocks=[convert_time(clock) for clock in Ω.f.clocks],
        values=[convert_rabi_amp(value) for value in Ω.f.values]
    )

    Δ = piecewise_linear(;
        clocks=[convert_time(clock) for clock in Δ.f.clocks],
        values=[convert_rabi_amp(value) for value in Δ.f.values]
    )
    return (ϕ,Ω,Δ,Δ_local)
end

parse_analog_rydberg_fields(ϕ,Ω,Δ) = throw(UndefVarError("Unable to parse Rydberg coefficients for Schema conversion, please use Real/Nothing for constant coefficients and Waveforms dynamic coefficients."))
# function will extract parameters, check of waveforms are compatible with IR
# then descretize waveforms and return them to be parsed into 
# effective Hamiltonian. 

function parse_analog_rydberg_params(h,params)
    atoms,ϕ,Ω,Δ = get_rydberg_params(h)
    # dispatch based on types because there are too many cases 
    # and each case requires a lot of code.
    ϕ,Ω,Δ,Δ_local = parse_analog_rydberg_fields(ϕ,Ω,Δ,params)

    return (atoms,ϕ,Ω,Δ,Δ_local)
end

