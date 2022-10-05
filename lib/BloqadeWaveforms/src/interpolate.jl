function piecewise_linear_interpolate(wf::Waveform{PiecewiseLinear{T,Interp},T}; 
    max_slope::Real=Inf64, 
    min_step::Real=0.0, 
    atol::Real = 1.0e-5) where {T<:Real,Interp}

    if atol < 0
        @warn "negative tolerance provided, taking absolute value."
        atol *= -1
    end

    if atol == 0 && ( max_slope == Inf64 || min_step == 0)
        error("Interpolation requires either a tolerance constraint or a slope and step constraint.")
    end

    c = wf.f.clocks
    v = wf.f.values

    @inbounds for i in 1:length(c)-1

        lb = c[i]
        ub = c[i+1]
        f_lb = v[i]
        f_ub = v[i+1]
        
        slope = (f_ub-f_lb)/(ub-lb)
        interval = ub-lb



        if abs(slope) > max_slope
            error("Waveform slope larger than constraint.")
        end

        if interval < min_step
            error("Waveform step smaller than constraint.")
        end
    end

    return wf

end

function piecewise_linear_interpolate(wf::Waveform{PiecewiseConstant{T},T}; 
    max_slope::Real=Inf64, 
    min_step::Real=0.0, 
    atol::Real = 1.0e-5) where {T<:Real}

    if atol < 0
        @warn "negative tolerance provided, taking absolute value."
        atol *= -1
    end

    if atol == 0 && ( max_slope == Inf64 || min_step == 0)
        error("Interpolation requires either a tolerance constraint or a slope and step constraint.")
    end

    c = wf.f.clocks
    v = wf.f.values
    clocks = Float64[c[1]]
    values = Float64[v[1]]


    itol = atol / (length(v) - 1) # distribute error over each step

    @inbounds for i in 1:length(v)-1
        t0 = c[i+1]        
        v[i+1] ≈ v[i] && continue # ignore steps which are small

        Δv = v[i+1]-v[i]

        Δt = if itol > 0
            # let error determine step size, error is equal to area between step and line functions
            4*itol/abs(Δv) 
        else
            # let max_slope / min_step determine step size to minimize error
            max(min_step,abs(Δv)/max_slope)
        end
        slope = Δv/Δt
        
        if abs(slope) > max_slope
            error("Requested tolerance for interpolation violates the slope constraint.")
        elseif Δt < min_step
            error("Requested tolerance for interpolation violates the step size constraint.")
        end

        if t0-Δt/2 < clocks[end]
            error("Distance between steps in piecewise constant waveform are too small to convert to piecewise linear.")
        end

        push!(values,v[i])
        push!(values,v[i+1])
        push!(clocks,t0-Δt/2)
        push!(clocks,t0+Δt/2)

    end

    push!(values,v[end])
    push!(clocks,c[end])


    return piecewise_linear(;clocks=clocks,values=values)

end

"""
    piecewise_linear_interpolate(waveform;[max_vale=Inf64,max_slope=Inf64,min_step=0.0])

Function which takes a waveform and translates it to a linear interpolation subject to some constraints. The function returns a piecewise linear waveform. if the Waveform is piecewise linear only the constraints will be checked. 

# Arguments

- `waveform`: ['Waveform'](@ref)  to be discretized.

# Keyword Arguments
- `min_step`: minimum possible step used in interpolation
- `max_slope`: Maximum possible slope used in interpolation
- `atol`: tolerance of interpolation, this is a bound to the area between the linear interpolation and the waveform.
"""
function piecewise_linear_interpolate(wf::Waveform; 
    max_slope::Real=Inf64, 
    min_step::Real=0.0, 
    atol::Real = 1.0e-5)
    
    if atol < 0
        @warn "negative tolerance provided, taking absolute value."
        atol *= -1
    end

    if atol == 0 && ( max_slope == Inf64 || min_step == 0)
        error("Interpolation requires either a tolerance constraint or a slope and step constraint.")
    end
    
    stack = NTuple{2,Float64}[(0.0,wf.duration)]
    intervals = NTuple{4,Float64}[]

    while !isempty(stack)
        lb,ub = pop!(stack)

        interval = (ub-lb)
        f_lb = wf(lb)
        f_ub = wf(ub)

        slope = (f_ub-f_lb)/interval

        lin_f = t -> slope .* (t .- lb) .+ f_lb

        area,_ = quadgk(t -> abs.(lin_f(t) .- wf.(t)),lb,ub,atol=eps())
        error_bound = max(atol*max(atol,interval),eps())

        if area*wf.duration > error_bound
            mid = (ub+lb)/2
            f_mid = wf(mid)

            next_slope = 2*max(abs(f_mid - f_lb),abs(f_ub - f_mid))/interval
            # only throw error if tolerance is non-zero
            # if tolerance is 0 simply stop adding to stack
            if next_slope > max_slope 
                atol > 0 && error("Requested tolerance for interpolation violates the slope constraint.")
                push!(intervals,(lb,ub,f_lb,slope))
            elseif interval < 2*min_step
                atol > 0 && error("Requested tolerance for interpolation violates the step size constraint.")
                push!(intervals,(lb,ub,f_lb,slope))
            else
                push!(stack,(mid,ub))
                push!(stack,(lb,mid))
            end

        else
            push!(intervals,(lb,ub,f_lb,slope))
        end

    end

    intervals = sort(intervals,by=ele->ele[1])

    clocks = Float64[each[1] for each in intervals]
    values = Float64[each[3] for each in intervals]

    push!(values,wf(wf.duration))
    push!(clocks,wf.duration)

    return piecewise_linear(;clocks=clocks,values=values)

end

function piecewise_constant_interpolate(wf::Waveform; 
    min_step::Real=0.0, 
    atol::Real = 1.0e-5)
    
    if atol < 0
        @warn "negative tolerance provided, taking absolute value."
        atol *= -1
    end

    if atol == 0 && min_step == 0
        error("Interpolation requires either a tolerance constraint or a minimum step constraint.")
    end
    
    
    area,_ = quadgk(t->wf.(t),0,wf.duration,atol=eps())
    value = area/wf.duration
    global_error,_ = quadgk(t -> abs.(value .- wf.(t)),0,wf.duration,atol=eps())

    chunks = Dict{NTuple{3,Float64},Float64}((global_error,0.0,wf.duration)=>value)


    while global_error > atol

        # split interval that can be split in half and has the largest error
        current_keys = collect(keys(chunks))
        sort!(current_keys,by=ele->ele[1])

        while length(current_keys) > 0
            _,lb,ub = current_keys[end]

            if (ub-lb) < 2*min_step 
                popfirst!(current_keys)
            else
                break
            end
        end

        if length(current_keys) == 0
            atol > 0 && error("Requested tolerance for interpolation violates the step size constraint.")
            break
        end

        max_key = current_keys[end]
        
        old_value = pop!(chunks,max_key)
        (old_error,lb,ub) = max_key
        mid = (ub+lb)/2
        total_area = old_value*(ub-lb)
        
        # only integrate over half of interval and use 
        # old area to calculate other part because lower_area + upper_area = total_area

        lower_area,_ = quadgk(wf,lb,mid,atol=eps())
        upper_area = total_area - lower_area

        lower_value = lower_area/(mid-lb)
        upper_value = upper_area/(ub-mid)

        lower_error,_ = quadgk(t -> abs.(lower_value .- wf.(t)),lb,mid,atol=eps())
        upper_error,_ = quadgk(t -> abs.(upper_value .- wf.(t)),mid,ub,atol=eps())

        chunks[(lower_error,lb,mid)] = lower_value
        chunks[(upper_error,mid,ub)] = upper_value

        # update global error
        global_error = (global_error - old_error) + lower_error + upper_error

    end

    
    current_keys = collect(keys(chunks))
    sort!(current_keys,by=ele->ele[2])

    clocks = Float64[]
    values = Float64[]
    for k in current_keys
        (error,lb,ub) = k

        push!(clocks,lb)
        push!(values,chunks[k])
    end

    push!(clocks,wf.duration)

    return piecewise_constant(;clocks=clocks,values=values)
end

function piecewise_constant_interpolate(wf::Waveform{PiecewiseConstant{T},T}; 
    min_step::Real=0.0, 
    atol::Real = 1.0e-5) where {T<:Real}

    if atol < 0
        @warn "negative tolerance provided, taking absolute value."
        atol *= -1
    end

    if atol == 0 && min_step == 0
        error("Interpolation requires either a tolerance constraint or a step constraint.")
    end

    c = wf.f.clocks

    @inbounds for i in 1:length(c)-1

        lb = c[i]
        ub = c[i+1]
        
        interval = ub-lb

        if interval < min_step
            error("Waveform step smaller than constraint.")
        end
    end

    return wf

end