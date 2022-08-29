

function get_rydberg_params(h::BloqadeExpr.RydbergHamiltonian)
    # extracts parameters from RydbergHamiltonian
    ϕ = nothing
    Ω = nothing
    Δ = nothing

    if h.rabi_term isa BloqadeExpr.SumOfX
        Ω = h.rabi_term.Ω.f # Ω is an instance of DivByTwo which has field f which is the function being divided by 2. 
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

    return (param,duration)
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

    # clock_samples = LinRange(0,duration,100)
    clock_samples = sort(duration .* rand(100))
    
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
    if any(s[2:end] .> s[1]*length(s)*eps(eltype(s)))
        throw(ErrorException("Local detuning waveforms cannot be decomposed into a product: Δ(i)⋅Δ(t)."))
    end 
    
    # use U vector to get the scal
    amplitude = (u * Diagonal(s))[:,1]
    i = argmax(abs.(amplitude))
    amplitude ./= amplitude[i]

    return (amplitude,param[i],duration)
    
end

# ϕ has a combination of both Ω and Δ constraints. 
# e.g. has to be global, but can begin and end on non-zero values.
function parse_static_rydberg_ϕ(param::Real,duration::Real,max_slope::Real,min_step::Real)
    _,ϕ = parse_static_rydberg_Δ(param,duration,max_slope,min_step)
    return ϕ
end

function parse_static_rydberg_ϕ(param::Nothing,duration::Real,max_slope::Real,min_step::Real)
    _,ϕ = parse_static_rydberg_Δ(param,duration,max_slope,min_step)
    return ϕ
end

function parse_dynamic_rydberg_ϕ(param::Waveform{F,T};duration=nothing) where {F,T<:Real}
    _,ϕ,duration = parse_dynamic_rydberg_Δ(param;duration)
    return (ϕ,duration)
end

parse_static_rydberg_ϕ(param::Vector{<:Real},duration::Real,max_slope::Real,min_step::Real) = parse_static_rydberg_Ω(param,duration,max_slope,min_step)
parse_dynamic_rydberg_ϕ(param::Vector{Waveform{F,T}};duration=nothing) where {F,T<:Real} = parse_dynamic_rydberg_Ω(param;duration)


@inline convert_units(value::Real,from,to) = uconvert(to,Quantity(value,from)).val
function convert_units(x::AbstractArray{S},from,to) where {S<:Real}
    y = similar(x)
    @inbounds for i in eachindex(x)
        y[i] = convert_units(x[i],from,to)
    end
    return y
end

function get_constraints(params)
    min_step = convert_units(params.rabi_time_min_step,s,μs)
    ϕ_max_slope = convert_units(params.rabi_frequency_phase_max_slope,rad/s,rad/μs)
    Ω_max_slope = convert_units(params.rabi_frequency_amplitude_max_slope,rad/s^2,rad*MHz/μs)
    Δ_max_slope = convert_units(params.rabi_detuning_max_slope,rad/s^2,rad*MHz/μs)
    return (min_step,ϕ_max_slope,Ω_max_slope,Δ_max_slope)
end

# no dynamic parameeters must throw error because `duration` can't be determined
function parse_analog_rydberg_fields(ϕ::ConstantParam,Ω::ConstantParam,Δ::ConstantParam,params) 
    throw(ErrorException("Schema requires at least one Waveform field."))
end

# one dynamic argument
function parse_analog_rydberg_fields(ϕ::DynamicParam,Ω::ConstantParam,Δ::ConstantParam,params) 
    (min_step,ϕ_max_slope,Ω_max_slope,Δ_max_slope) = get_constraints(params)

    ϕ,duration = parse_dynamic_rydberg_ϕ(ϕ)
    Ω = parse_static_rydberg_Ω(Ω,duration,Ω_max_slope,min_step)
    Δ_local,Δ = parse_static_rydberg_Δ(Δ,duration,Δ_max_slope,min_step)

    ϕ = BloqadeWaveforms.discretize(ϕ;
        max_slope=ϕ_max_slope,
        min_step=min_step,
        tol=params.waveform_tolerance,
    )

    ϕ = piecewise_linear(;
        clocks=convert_units(ϕ.f.clocks,μs,s),
        values=convert_units(ϕ.f.values,rad,rad)
    )

    Ω = piecewise_linear(;
        clocks=convert_units(Ω.f.clocks,μs,s),
        values=convert_units(Ω.f.values,rad*MHz,rad/s)
    )

    Δ = piecewise_linear(;
        clocks=convert_units(Δ.f.clocks,μs,s),
        values=convert_units(Δ.f.values,rad*MHz,rad/s)
    )
    return (ϕ,Ω,Δ,Δ_local)
end

function parse_analog_rydberg_fields(ϕ::ConstantParam,Ω::DynamicParam,Δ::ConstantParam,params)
    (min_step,ϕ_max_slope,Ω_max_slope,Δ_max_slope) = get_constraints(params)

    Ω,duration = parse_dynamic_rydberg_ϕ(Ω)
    ϕ = parse_static_rydberg_ϕ(ϕ,duration,ϕ_max_slope,min_step)
    Δ_local,Δ = parse_static_rydberg_Δ(Δ,duration,Δ_max_slope,min_step)

    Ω = BloqadeWaveforms.discretize(Ω;
        max_slope=Ω_max_slope,
        min_step=min_step,
        tol=params.waveform_tolerance,
    )

    ϕ = piecewise_linear(;
        clocks=convert_units(ϕ.f.clocks,μs,s),
        values=convert_units(ϕ.f.values,rad,rad)
    )

    Ω = piecewise_linear(;
        clocks=convert_units(Ω.f.clocks,μs,s),
        values=convert_units(Ω.f.values,rad*MHz,rad/s)
    )

    Δ = piecewise_linear(;
        clocks=convert_units(Δ.f.clocks,μs,s),
        values=convert_units(Δ.f.values,rad*MHz,rad/s)
    )
    return (ϕ,Ω,Δ,Δ_local)
end

function parse_analog_rydberg_fields(ϕ::ConstantParam,Ω::ConstantParam,Δ::DynamicParam,params) 
    (min_step,ϕ_max_slope,Ω_max_slope,Δ_max_slope) = get_constraints(params)

    Δ_local,Δ,duration = parse_dynamic_rydberg_Δ(Δ)
    ϕ = parse_static_rydberg_ϕ(ϕ,duration,ϕ_max_slope,min_step)
    Ω = parse_static_rydberg_Ω(Ω,duration,Ω_max_slope,min_step)

    Δ = BloqadeWaveforms.discretize(Δ;
        max_slope=Δ_max_slope,
        min_step=min_step,
        tol=params.waveform_tolerance,
    )

    ϕ = piecewise_linear(;
        clocks=convert_units(ϕ.f.clocks,μs,s),
        values=convert_units(ϕ.f.values,rad,rad)
    )

    Ω = piecewise_linear(;
        clocks=convert_units(Ω.f.clocks,μs,s),
        values=convert_units(Ω.f.values,rad*MHz,rad/s)
    )

    Δ = piecewise_linear(;
        clocks=convert_units(Δ.f.clocks,μs,s),
        values=convert_units(Δ.f.values,rad*MHz,rad/s)
    )
    return (ϕ,Ω,Δ,Δ_local)
end

# two dynamic arguments
function parse_analog_rydberg_fields(ϕ::DynamicParam,Ω::DynamicParam,Δ::ConstantParam,params) 
    (min_step,ϕ_max_slope,Ω_max_slope,Δ_max_slope) = get_constraints(params)

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
        clocks=convert_units(ϕ.f.clocks,μs,s),
        values=convert_units(ϕ.f.values,rad,rad)
    )

    Ω = piecewise_linear(;
        clocks=convert_units(Ω.f.clocks,μs,s),
        values=convert_units(Ω.f.values,rad*MHz,rad/s)
    )

    Δ = piecewise_linear(;
        clocks=convert_units(Δ.f.clocks,μs,s),
        values=convert_units(Δ.f.values,rad*MHz,rad/s)
    )
    return (ϕ,Ω,Δ,Δ_local)
end

function parse_analog_rydberg_fields(ϕ::DynamicParam,Ω::ConstantParam,Δ::DynamicParam,params) 
    (min_step,ϕ_max_slope,Ω_max_slope,Δ_max_slope) = get_constraints(params)

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
        clocks=convert_units(ϕ.f.clocks,μs,s),
        values=convert_units(ϕ.f.values,rad,rad)
    )

    Ω = piecewise_linear(;
        clocks=convert_units(Ω.f.clocks,μs,s),
        values=convert_units(Ω.f.values,rad*MHz,rad/s)
    )

    Δ = piecewise_linear(;
        clocks=convert_units(Δ.f.clocks,μs,s),
        values=convert_units(Δ.f.values,rad*MHz,rad/s)
    )
    return (ϕ,Ω,Δ,Δ_local)
end

function parse_analog_rydberg_fields(ϕ::ConstantParam,Ω::DynamicParam,Δ::DynamicParam,params) 
    (min_step,ϕ_max_slope,Ω_max_slope,Δ_max_slope) = get_constraints(params)

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
        clocks=convert_units(ϕ.f.clocks,μs,s),
        values=convert_units(ϕ.f.values,rad,rad)
    )

    Ω = piecewise_linear(;
        clocks=convert_units(Ω.f.clocks,μs,s),
        values=convert_units(Ω.f.values,rad*MHz,rad/s)
    )

    Δ = piecewise_linear(;
        clocks=convert_units(Δ.f.clocks,μs,s),
        values=convert_units(Δ.f.values,rad*MHz,rad/s)
    )
    return (ϕ,Ω,Δ,Δ_local)
end

# three dynamic arguments
function parse_analog_rydberg_fields(ϕ::DynamicParam,Ω::DynamicParam,Δ::DynamicParam,params) 
    (min_step,ϕ_max_slope,Ω_max_slope,Δ_max_slope) = get_constraints(params)

    ϕ,duration = parse_dynamic_rydberg_ϕ(ϕ)
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
        clocks=convert_units(ϕ.f.clocks,μs,s),
        values=convert_units(ϕ.f.values,rad,rad)
    )

    Ω = piecewise_linear(;
        clocks=convert_units(Ω.f.clocks,μs,s),
        values=convert_units(Ω.f.values,rad*MHz,rad/s)
    )

    Δ = piecewise_linear(;
        clocks=convert_units(Δ.f.clocks,μs,s),
        values=convert_units(Δ.f.values,rad*MHz,rad/s)
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

    atoms = map(atoms) do pos 
                return (convert_units(pos[1],μm,m),convert_units(pos[2],μm,m))
            end
    ϕ,Ω,Δ,Δ_local = parse_analog_rydberg_fields(ϕ,Ω,Δ,params)

    return (atoms,ϕ,Ω,Δ,Δ_local)
end

