# this function is required to calculate difference
# between waveforms which have different durations.
# the behavior is to pad the shorter waveform at the
# end with 0. This is useful when the user passes
# in a duration that doesn't neatly fit into the
# time resolution of the device. In which case 
# this is reqiured to check the error.
function norm_diff_durations(A,B)

    (A,B) = ( A.duration > B.duration ? (A,B) : (B,A))
    A = if A.duration > B.duration
        T_diff = A.duration - B.duration
        append(A,constant(;duration=T_diff,value=0))
    end
    return norm(A-B)
end

# simple function, round number to the nearest multiple of `res`
function set_resolution(x::Real,res::Real)
    @assert res > 0
    return round(Int(round(x / res)) * res,sigdigits=15)
end

# throws error for none-scalar fields
function  check_global(field,name) 
    field isa Vector && error("Failed to transform $name to hardware, $name must be global drive.")
end

# checks argument to make sure it is a waveform or a container of waveforms.
function check_waveform(field,name)
    field isa Waveform && return 

    (isnothing(field) || field isa Number || is_time_function(field)) && error("Failed to transform $name to hardware, $name must be a Waveform.\n If $name is constant or zero use `BloqadeWaveforms.constant`.")
    field isa Vector{<:Number} && error("Failed to transform $name to hardware, $name must be a vector Waveforms.\n If $name contains constant values use `BloqadeWaveforms.constant`.")
    if field isa Vector 
        any( !(ele isa Waveform) for ele in field) && error("Failed to transform $name to hardware, $name must be a Waveform.")

        duration = field[1].duration

        any(duration!=ele.duration for ele in field) && error("Failed to transform $name to hardware, all Waveforms must have the same duration.")
        
    end
end
# warns user if the duration of the waveform will be rounded by the time resolution.
function warn_duration(time_res,field,name)
    duration = if field isa Vector
        field[1].duration
    else
        field.duration
    end
    if !(set_resolution(duration,time_res) ≈ duration)
        @info "waveform $name duration will be rounded during hardware transformation."
    end
end

# TODO: move this to BloqadeWaveforms
# function inserts ramps at the beginning and end to set correct initial and final values
function pin_waveform_edges(wf::Waveform,
    max_slope::Real,
    begin_value::Real,
    end_value::Real)

    duration = wf.duration



    t_begin = if !isapprox(Ω(0.0), begin_value;atol=eps(),rtol=√eps())
        ramp_up =  (sign(Ω(0.0)-begin_value)*max_slope)
        lin_ramp_begin = Waveform(t -> ramp_up .* t .+ begin_value, duration)
        find_zero(wf-lin_ramp_begin,(0.0,duration))
    else
        0.0
    end

    t_end = if !isapprox(Ω(duration), end_value;atol=eps(),rtol=√eps())
        ramp_down = (sign(Ω(duration)-end_value)*max_slope)
        lin_ramp_end = Waveform(t -> ramp_down .* (t.-duration) .+ end_value, duration)
        find_zero(wf-lin_ramp_end,(0.0,duration))
    else
        duration
    end
        
    if t_begin > 0 && t_end < duration
        mid_wf = Waveform(t->wf.f(t.+t_begin),duration - t_begin - t_end)

        start_wf = linear_ramp(;duration=t_begin,
            start_value=begin_value,
            stop_value=wf(t_begin)
        )

        end_wf = linear_ramp(;duration=duration-t_end,
            start_value=wf(t_end),
            end_value=end_value
        )

        return append(start_wf,mid_wf,end_wf)
    elseif t_begin > 0 
        end_wf = Waveform(t->wf.f(t.+t_begin),duration - t_begin)

        start_wf = linear_ramp(;duration=t_begin,
            start_value=begin_value,
            stop_value=wf(t_begin)
        )

        return append(start_wf,end_wf)
    elseif t_end < duration
        start_wf = Waveform(wf.f,duration - t_end)

        end_wf = linear_ramp(;duration=duration-t_end,
            start_value=wf(t_end),
            end_value=end_value
        )

        return append(start_wf,end_wf)
    else
        return wf      
    end

    return new_wf

end
# TODO: move this to BloqadeWaveforms
# does the SVD method to decompose local drives into an outer produce of masks and functions.
function find_local_masks(values::Array{T,2},name;ntrunc=1,assert_truncation=false) where T
    
    (l,w) = size(values)
    
    # project out uniform amplitude waveforms
    ones_norm = ones(w,1)./sqrt(w)
    avg = (values*ones_norm)*transpose(ones_norm)
    results = []
    values = values - avg
    avg = avg[:,1]

    if ntrunc > 0
        # use SVD to decompose the non-uniform waveforms
        # if there are more than one singular value larger than 0 the remaining pattern
        # can't be described as a local detuning mask times a single time-dependent function.
        # here we truncate the waveform only keeping the first singular value.
        u,s,v_T = svd(values)
        contained_weight = sum(s[:ntrunc])
        remaining_weight = sum(s[ntrunc+1:end])

        if remaining_weight > eps()*length(s)*contained_weight
            error_or_warn(assert_truncation,"Cannot decompose waveform $name into product of masks and scalar functions.")
        end

        for i in 1:ntrunc

            if all(u[:,i] .< 0) # make u all positive
                u[:,i] .*= -1
                mask = - s[i] .* v_T[:,i]
            end

            min = minimum(mask)
            max = maximum(mask)

            mask = (mask .- min) ./ (max - min)

            f = u[:,i] .* (max - min)
            avg .+= u[:,i] .* min

            push!(results,(f,mask))
        end
    end

    push!(results,(avg,ones(w)))
    return results

end

function hardware_transform_Ω(Ω,device_capabilities::DeviceCapabilities)
    time_res = device_capabilities.rydberg.global_value.timeResolution
    min_step = device_capabilities.rydberg.global_value.timeDeltaMin
    rabi_res = device_capabilities.rydberg.global_value.rabiFrequencyResolution
    max_slope = device_capabilities.rydberg.global_value.rabiFrequencySlewRateMax

    check_waveform(Ω,:Ω)
    check_global(Ω,:Ω)
    warn_duration(time_res,Ω,:Ω)

    Ω = if !isapprox(Ω(0.0), 0.0;atol=eps(),rtol=√eps()) || !isapprox(Ω(duration),0.0;atol=eps(),rtol=√eps())
        @info "During hardware transform: Ω start and/or end values are not 0. adding ramp(s) to fix endpoints."
        pin_waveform_edges(Ω,max_slope,0,0)
    end

    Ωt = if Ω isa PiecewiseLinearWaveform
        piecewise_linear(;
            clocks=set_resolution.(Ω.f.clocks,time_res),
            values=set_resolution.(Ω.f.values,rabi_res)
        )

    elseif Ω isa Waveform
        
        Ω_interp = piecewise_linear_interpolate(Ω,max_slope=max_slope,min_step=min_step)
        piecewise_linear(;
            clocks=set_resolution.(Ω_interp.f.clocks,time_res),
            values=set_resolution.(Ω_interp.f.values,rabi_res)
        )
    end

    return Ωt,norm_diff_durations(Ω,Ωt)
end

function hardware_transform_ϕ(ϕ,device_capabilities::DeviceCapabilities)
    time_res = device_capabilities.rydberg.global_value.timeResolution
    min_step = device_capabilities.rydberg.global_value.timeDeltaMin
    phase_res = device_capabilities.rydberg.global_value.phaseResolution
    max_slope = device_capabilities.rydberg.global_value.phaseSlewRateMax
    check_global(ϕ,:ϕ)
    check_waveform(ϕ,:ϕ)
    warn_duration(time_res,ϕ,:ϕ)

    ϕt = if ϕ isa PiecewiseLinearWaveform
        piecewise_linear(;
            clocks=set_resolution.(ϕ.f.clocks,time_res),
            values=set_resolution.(ϕ.f.values,phase_res)
        )        
    elseif ϕ isa Waveform
        # arbitrary waveform must transform
        
        ϕ_interp = piecewise_linear_interpolate(ϕ,max_slope=max_slope,min_step=min_step)
        piecewise_linear(;
            clocks=set_resolution.(ϕ_interp.f.clocks,time_res),
            values=set_resolution.(ϕ_interp.f.values,phase_res)
        )
    end

    return ϕt,norm_diff_durations(ϕ,ϕt)
end

function hardware_transform_Δ(Δ,device_capabilities::DeviceCapabilities)

    time_res = device_capabilities.rydberg.global_value.timeResolution
    min_step = device_capabilities.rydberg.global_value.timeDeltaMin
    
    detune_res = device_capabilities.rydberg.global_value.phaseResolution
    local_res = device_capabilities.rydberg.local_value.localDetuningResolution
    common_detune_res = device_capabilities.rydberg.local_value.commonDetuningResolution
    max_slope = device_capabilities.rydberg.global_value.phaseSlewRateMax

    check_waveform(Δ,:Δ)
    warn_duration(time_res,Δ,:Δ)


    if Δ isa PiecewiseLinearWaveform
        Δt = piecewise_linear(;
            clocks=set_resolution.(Δ.f.clocks,time_res),
            values=set_resolution.(Δ.f.values,detune_res)
        )
        Δ_mask = (Δ=Δt,δ=nothing,Δi=1.0)
        
        return Δt,norm_diff_durations(Δ,Δt),Δ_mask

    elseif Δ isa Waveform
        # arbitrary waveform must transform
        
        Δ_interp = piecewise_linear_interpolate(Δ,max_slope=max_slope,min_step=min_step)
        Δt = piecewise_linear(;
            clocks=set_resolution.(Δ_interp.f.clocks,time_res),
            values=set_resolution.(Δ_interp.f.values,detune_res)
        )
        Δ_mask = (Δ=Δt,δ=nothing,Δi=1.0)

        return Δt,norm_diff_durations(Δ,Δt),Δ_mask
    elseif Δ isa Vector
        
        nsteps = Int(Δ[1].duration÷min_step)
        clocks = collect(LinRange(0.0,duration,nsteps))
        clocks = set_resolution.(clocks,time_res)

        values = zeros(length(clocks),length(Δ))

        @inbounds for (j,δ) in enumerate(Δ)
            values[:,j] = δ(clocks)
        end

        ((δ_values,Δi),(Δ_values,_)) = find_local_masks(values,:Δ)

        Δ_values = set_resolution.(Δ_values,detune_res)
        δ_values = set_resolution.(δ_values,common_detune_res)
        Δi = set_resolution.(Δi,local_res)

        Δ_mask = (
            Δ=piecewise_linear(;clocks=clocks,values=Δ_values),
            δ=piecewise_linear(;clocks=clocks,values=Δ_values),
            Δi=Δi
        )

        error = [norm_diff_durations(Δ[i],Δt+δi*δ) for (i,δi) in enumerate(Δi)]
        Δt = [piecewise_linear(;clocks=clocks,values=Δ_values+δi.*δ_values) for δi in Δi]

        return Δt,error,Δ_mask
    end


end

function hardware_transform_parse(H::BloqadeExpr.RydbergHamiltonian;device_capabilities::DeviceCapabilities=get_device_capabilities())
    (atoms,ϕ,Ω,Δ) = get_rydberg_params(H)

    ϕ,ϕ_error = hardware_transform_ϕ(ϕ,device_capabilities)
    Ω,Ω_error = hardware_transform_Ω(Ω,device_capabilities)
    Δ,Δ_error,Δ_mask = hardware_transform_Δ(Δ,device_capabilities)
    info = (ϕ=ϕ_error,Ω=Ω_error,Δ=Δ_error,Δ_mask=Δ_mask)

    return (ϕ,Ω,Δ,info)
end

# public API exposed here: 
function hardware_transform(H::BloqadeExpr.RydbergHamiltonian;device_capabilities::DeviceCapabilities=get_device_capabilities())
    ϕ,Ω,Δ,info = hardware_transform_parse(H;device_capabilities=device_capabilities)
    hardware_H = rydberg_h(atoms,ϕ=ϕ,Ω=Ω,Δ=Δ)

    return hardware_H,info
end
