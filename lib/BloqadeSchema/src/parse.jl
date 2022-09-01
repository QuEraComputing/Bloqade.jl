

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
    max_time = convert_units(params.rabi_time_maximum_value,s,μs)
    ϕ_max_slope = convert_units(params.rabi_frequency_phase_max_slope,rad/s,rad/μs)
    Ω_max_slope = convert_units(params.rabi_frequency_amplitude_max_slope,rad/s^2,rad*MHz/μs)
    Δ_max_slope = convert_units(params.rabi_detuning_max_slope,rad/s^2,rad*MHz/μs)
    return (min_step,max_time,ϕ_max_slope,Ω_max_slope,Δ_max_slope)
end

# parsing individual fields
# we split them up because each 
# field has different constraints. 

function parse_static_rydberg_Ω(Ω::ConstantParam,duration::Real,params)
    if isnothing(Ω) || Ω == 0
        return piecewise_linear(;clocks=Float64[0.0,duration],values=Float64[0.0,0.0])
    elseif Ω isa Real
        constraints = get_constraints(params)
        min_step = constraints[1]
        max_slope = constraints[2]

        step = max(min_step,abs(Ω/max_slope))

        clocks = Float64[0.0,step,duration-step,duration]
        values = Float64[0.0,Ω,Ω,0.0]
        waveform = piecewise_linear(;clocks,values)
    
        return waveform
    else
        throw(ErrorException("Rabi field amplitude must be global drive."))
    end
end

function parse_static_rydberg_Δ(Δ::ConstantParam,duration::Real)
    if isnothing(Δ)
        return (1.0,nothing,piecewise_linear(;clocks=Float64[0.0,duration],values=Float64[0.0,0.0]))
    elseif Δ isa Real
       return (1.0,nothing,piecewise_linear(;clocks=Float64[0.0,duration],values=Float64[Δ,Δ]))
    else
        Δi = Float64[val for val in Δ]
        Δ0 = minimum(Δi)
        Δi .-= Δ0

        δ_max = maximum(Δi)
        Δi ./= δ_max

        Δ = piecewise_linear(;clocks=Float64[0.0,duration],values=Float64[Δ0,Δ0])
        δ = piecewise_linear(;clocks=Float64[0.0,duration],values=Float64[δ_max,δ_max])

        return (Δi,δ,Δ)
    end
end

function parse_static_rydberg_ϕ(ϕ::ConstantParam,duration::Real)
    if isnothing(ϕ)
        return piecewise_linear(;clocks=Float64[0.0,duration],values=Float64[0.0,0.0])
    elseif ϕ isa Real
       return piecewise_linear(;clocks=Float64[0.0,duration],values=Float64[ϕ,ϕ])
    else
        throw(ErrorException("Rabi field phase must be global drive.")) 
    end
end

function parse_dynamic_rydberg_Ω(Ω::DynamicParam,params;duration=nothing)
    if Ω isa Waveform
        min_step,max_time,ϕ_max_slope,Ω_max_slope,Δ_max_slope = get_constraints(params)

        duration = (duration ≡ nothing ? Ω.duration : duration)
    
        if !(duration ≡ nothing) && Ω.duration != duration
            throw(ErrorException("Waveform durations do not match."))
        end
    
        if Ω.duration < min_step
            throw(ErrorException("Waveform Ω duration must be larger than minimum step, $(min_step) μs."))
        end

        if Ω.duration > max_time
            throw(ErrorException("Waveform Ω duration must be shorter than maximum time, $(max_time) μs."))
        end

        if !isapprox(Ω(0.0), 0.0;atol=eps(),rtol=eps()) || !isapprox(Ω(duration),0.0;atol=eps(),rtol=eps())
            throw(ErrorException("Rabi Drive must start and end with value 0."))
        end

        Ω = BloqadeWaveforms.discretize(Ω;
            max_slope=Ω_max_slope,
            min_step=min_step,
            tol=params.waveform_tolerance,
        )

        return Ω,duration
    else
        throw(ErrorException("Rabi field amplitude must be global drive."))
    end
end

function parse_dynamic_rydberg_ϕ(ϕ::DynamicParam,params;duration=nothing)
    if ϕ isa Waveform
        min_step,max_time,ϕ_max_slope,Ω_max_slope,Δ_max_slope = get_constraints(params)

        duration = (duration ≡ nothing ? ϕ.duration : duration)
    
        if !(duration ≡ nothing) && ϕ.duration != duration
            throw(ErrorException("Waveform durations do not match."))
        end
    
        if ϕ.duration < min_step
            throw(ErrorException("Waveform ϕ duration must be larger than minimum step, $(min_step) μs."))
        end

        if ϕ.duration > max_time
            throw(ErrorException("Waveform ϕ duration must be shorter than maximum time, $(max_time) μs."))
        end

        ϕ = BloqadeWaveforms.discretize(ϕ;
            max_slope=ϕ_max_slope,
            min_step=min_step,
            tol=params.waveform_tolerance,
        )

        return ϕ,duration
    else
        throw(ErrorException("Rabi field phase must be global drive."))
    end
        
end


function parse_dynamic_rydberg_Δ(Δ::DynamicParam,params;duration=nothing)
    min_step,max_time,ϕ_max_slope,Ω_max_slope,Δ_max_slope = get_constraints(params)

    if Δ isa Waveform
        duration = (duration ≡ nothing ? Δ.duration : duration)
    
        if !(duration ≡ nothing) && Δ.duration != duration
            throw(ErrorException("Waveform durations do not match."))
        end

        if Δ.duration < min_step
            throw(ErrorException("Waveform Δ duration must be larger than minimum step, $(min_step) μs."))
        end

        if Δ.duration > max_time
            throw(ErrorException("Waveform Δ duration must be shorter than maximum time, $(max_time) μs."))
        end

        Δ = BloqadeWaveforms.discretize(Δ;
            max_slope=Δ_max_slope,
            min_step=min_step,
            tol=params.waveform_tolerance,
        )

        return (1.0,nothing,Δ,duration)
    else
        durations = [f.duration for f in Δ]
    
        duration = (duration ≡ nothing ? durations[1] : duration)
    
        if !(duration ≡ nothing) && !all(duration .== durations)
            throw(ErrorException("Waveform durations do not match."))
        end
      
        if durations[1] < min_step
            throw(ErrorException("Waveform Δ duration must be larger than minimum step, $(min_step) μs."))
        end

        if durations[1] > max_time
            throw(ErrorException("Waveform Δ duration must be shorter than maximum time, $(max_time) μs."))
        end


        nsteps = Int(durations[1]÷min_step)
        clocks = collect(LinRange(0.0,duration,nsteps))
        Δti = zeros(length(clocks),length(Δ))
    
        for (i,clock) in enumerate(clocks)
            for (j,f) in enumerate(Δ)
                Δti[i,j] = f(clock)
            end
        end
        
        # project out uniform amplitude waveforms
        ones_norm = ones(length(Δ),1)./sqrt(length(Δ))
        Δt = (Δti*ones_norm)*transpose(ones_norm)
        δti = Δti - Δt

        # use SVD to decompose the non-uniform waveforms
        # if there are more than one singular value larger than 0 the remaining pattern
        # can't be described as a local detuning mask times a single time-dependent function.
        u,s,v_T = svd(δti)
    
        if any(s[2:end] .> s[1]*length(s)*eps())
            throw(ErrorException("Local detuning waveforms cannot be decomposed into a product: Δ(i,t) = Δ(t)+ Δ_i⋅δ(t)."))
        end 
    
        # now take one column since they are all the same.
        Δt = Δt[:,1]
        # get the function values by multiplying s into u.
        # because there is only one non-zero singular value we 
        # just take the first column of u. 
        δt = s[1] .* u[:,1]
        # there are the local lattice site amplitudes.
        amps = v_T[:,1]
        
        # the lattice amplitudes must be within 0 and 1.
        amps_min = minimum(amps)
        amps .-= amps_min
        amps_max = maximum(amps)
        Δi = amps ./ amps_max

        # shift the minimum amplitude into the global drive.
        Δt .+= δt .* amps_min
    
        # rescale amplitudes to be a maximum of 1, absorb that scaling into δt
        δt .*= amps_max

        # get waveforms
        Δ = piecewise_linear(;
            clocks = clocks,
            values = Δt
        )

        δ = piecewise_linear(;
            clocks = clocks,
            values = δt
        )

        return (Δi,δ,Δ,duration)
    end
end


# no dynamic parameeters must throw error because `duration` can't be determined
function parse_analog_rydberg_fields(ϕ::ConstantParam,Ω::ConstantParam,Δ::ConstantParam,params) 
    throw(ErrorException("Schema requires at least one Waveform field."))
end

# one dynamic argument
function parse_analog_rydberg_fields(ϕ::DynamicParam,Ω::ConstantParam,Δ::ConstantParam,params) 

    ϕ,duration = parse_dynamic_rydberg_ϕ(ϕ,params)
    Ω = parse_static_rydberg_Ω(Ω,duration,params)
    Δ_local,δ,Δ = parse_static_rydberg_Δ(Δ,duration)


    return (ϕ,Ω,Δ,δ,Δ_local)
end

function parse_analog_rydberg_fields(ϕ::ConstantParam,Ω::DynamicParam,Δ::ConstantParam,params)

    Ω,duration = parse_dynamic_rydberg_Ω(Ω,params)
    ϕ = parse_static_rydberg_ϕ(ϕ,duration)
    Δ_local,δ,Δ = parse_static_rydberg_Δ(Δ,duration)

    return (ϕ,Ω,Δ,δ,Δ_local)
end

function parse_analog_rydberg_fields(ϕ::ConstantParam,Ω::ConstantParam,Δ::DynamicParam,params) 

    Δ_local,δ,Δ,duration = parse_dynamic_rydberg_Δ(Δ,params)
    ϕ = parse_static_rydberg_ϕ(ϕ,duration)
    Ω = parse_static_rydberg_Ω(Ω,duration,params)

    return (ϕ,Ω,Δ,δ,Δ_local)
end

# two dynamic arguments
function parse_analog_rydberg_fields(ϕ::DynamicParam,Ω::DynamicParam,Δ::ConstantParam,params) 

    ϕ,duration = parse_dynamic_rydberg_ϕ(ϕ,params)
    Ω,duration = parse_dynamic_rydberg_Ω(Ω,params;duration)
    Δ_local,δ,Δ = parse_static_rydberg_Δ(Δ,duration)

    return (ϕ,Ω,Δ,δ,Δ_local)
end

function parse_analog_rydberg_fields(ϕ::DynamicParam,Ω::ConstantParam,Δ::DynamicParam,params) 

    ϕ,duration = parse_dynamic_rydberg_ϕ(ϕ,params)
    Δ_local,δ,Δ,duration = parse_dynamic_rydberg_Δ(Δ,params;duration)
    Ω = parse_static_rydberg_Ω(Ω,duration,params)

    return (ϕ,Ω,Δ,δ,Δ_local)
end

function parse_analog_rydberg_fields(ϕ::ConstantParam,Ω::DynamicParam,Δ::DynamicParam,params) 

    Ω,duration = parse_dynamic_rydberg_Ω(Ω,params)
    Δ_local,δ,Δ,duration = parse_dynamic_rydberg_Δ(Δ,params;duration)
    ϕ = parse_static_rydberg_ϕ(ϕ,duration)

    return (ϕ,Ω,Δ,δ,Δ_local)
end

# three dynamic arguments
function parse_analog_rydberg_fields(ϕ::DynamicParam,Ω::DynamicParam,Δ::DynamicParam,params) 

    ϕ,duration = parse_dynamic_rydberg_ϕ(ϕ,params)
    Ω,duration = parse_dynamic_rydberg_Ω(Ω,params;duration)
    Δ_local,δ,Δ,duration = parse_dynamic_rydberg_Δ(Δ,params;duration)
    
    return (ϕ,Ω,Δ,δ,Δ_local)
end

parse_analog_rydberg_fields(ϕ,Ω,Δ) = throw(UndefVarError("Unable to parse Rydberg coefficients for Schema conversion, please use Real/Nothing for constant coefficients and Waveforms dynamic coefficients."))
# function will extract parameters, check of waveforms are compatible with IR
# then descretize waveforms and return them to be parsed into 
# effective Hamiltonian. 

function parse_analog_rydberg_params(h::BloqadeExpr.RydbergHamiltonian,params::SchemaConversionParams)
    atoms,ϕ,Ω,Δ = get_rydberg_params(h)
    # dispatch based on types because there are too many cases 
    # and each case requires a lot of code.

    atoms = map(atoms) do pos 
                return (convert_units(pos[1],μm,m),convert_units(pos[2],μm,m))
            end
    ϕ,Ω,Δ,δ,Δ_local = parse_analog_rydberg_fields(ϕ,Ω,Δ,params)

    return (atoms,ϕ,Ω,Δ,δ,Δ_local)
end

