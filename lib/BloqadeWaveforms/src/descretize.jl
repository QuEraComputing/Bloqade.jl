



function descretize_waveform(wf::Waveform{PiecewiseLinear{T,Interp},T}; 
    max_value::Real=Inf, 
    max_slope::Real=Inf, 
    min_step::Real=0, 
    tol::Float64 = 1.0e-3) where {T<:Real,Interp}

    c = wf.f.clocks
    v = wf.f.values

    @inbounds for i in 1:length(c)-1

        lb = c[i]
        ub = c[i+1]
        f_lb = v[i]
        f_ub = v[i+1]
        
        slope = (f_ub-f_lb)/(ub-lb)
        interval = ub-lb


        if abs(f_lb) > max_value || abs(f_ub) > max_value
            throw(ErrorException("Waveform exceeds maximum value."))
        end

        if interval < min_step
            throw(ErrorException("Waveform step smaller than constraint."))
        end

        if abs(slope) > max_slope
            throw(ErrorException("Waveform slope larger than constraint."))
        end
        
    end

    return wf

end

function descretize_waveform(wf::Waveform{PiecewiseConstant{T},T}; 
    max_value::Real=Inf, 
    max_slope::Real=Inf, 
    min_step::Real=0, 
    tol::Float64 = 1.0e-3) where {T<:Real}

    c = wf.f.clocks
    v = wf.f.values
    clocks = Float64[c[1]]
    values = Float64[v[1]]

    if any(abs(ele)>max_value for ele in v)
        throw(ErrorException("Waveform exceeds maximum value."))
    end

    itol = tol / (length(v) - 1) # distribute error over each step

    for i in 1:length(v)-1
        t0 = c[i+1]        
        Δv = v[i+1]-v[i]
        Δt = 4*itol/abs(Δv) # let error determine step size, error is equal to area between step and line functions
        slope = Δv/Δt
        
        if Δt < min_step || abs(slope) > max_slope
            throw(ErrorException("Descretization cannot obtain requested tolerance given the step size and slope constraints."))
        end

        if t0-Δt/2 < clocks[end]
            throw(ErrorException("Distance between steps in waveform are too small for requested descretization tolerance."))
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

function descretize_waveform(wf::Waveform; 
    max_value::Real=Inf, 
    max_slope::Real=Inf, 
    min_step::Real=0, 
    tol::Float64 = 1.0e-5)
    
    # TODO:
    # need to optimize stack and intervals

    stack = Tuple{Float64,Float64}[(0.0,wf.duration)]
    intervals = Tuple{Float64,Float64,Float64,Float64}[]
    wf_wrapper = t -> sample_values(wf,t)
    while !isempty(stack)
        lb,ub = pop!(stack)

        interval = (ub-lb)
        f_lb = wf(lb)
        f_ub = wf(ub)

        if abs(f_lb) > max_value || abs(f_ub) > max_value
            throw(ErrorException("Waveform exceeds maximum value."))
        end

        slope = (f_ub-f_lb)/interval

        lin_f = t -> slope .* (t .- lb) .+ f_lb

        area,_ = quadgk(t -> abs.(lin_f(t) .- wf_wrapper(t)),lb,ub)
        error_bound = tol*max(tol,interval)

        if area*wf.duration > error_bound
            mid = (ub+lb)/2
            f_mid = wf(mid)

            next_slope = 2*max(abs(f_mid - f_lb),abs(f_ub - f_mid))/interval
            
            if interval < 2*min_step || next_slope > max_slope
                # can't satisfy constraint. stop refining grid
                throw(ErrorException("Descretization cannot obtain requested tolerance given the step size and slope constraints."))
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