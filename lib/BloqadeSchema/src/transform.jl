# this function is required to calculate difference
# between waveforms which have different durations.
# the behavior is to pad the shorter waveform at the
# end with 0. This is useful when the user passes
# in a duration that doesn't neatly fit into the
# time resolution of the device. In which case 
# this is reqiured to check the error.
function norm_diff_durations(A::Waveform,B::Waveform)
    A.duration==B.duration && return norm(A-B)

    if A.duration > B.duration
        T_diff = A.duration - B.duration
        B = append(B,constant(;duration=T_diff,value=0))
        return norm(A-B)
    else
        T_diff = B.duration - B.duration
        A = append(A,constant(;duration=T_diff,value=0))
        return norm(A-B)
    end
end

# simple function, round number to the nearest multiple of `res`
function set_resolution(x::Real,res::Real)
    @assert res > 0
    return round(Int(round(x / res)) * res,sigdigits=14)
end

# throws error for none-scalar fields
function  check_global(field,name) 
    field isa AbstractArray && error("Failed to transform $name to hardware, $name must be global drive.")
end

# checks argument to make sure it is a waveform or a container of waveforms.
function check_waveform(field,name)
    field isa Waveform{F,T} where {F,T<:Real} && return 
        
    field isa Vector{<:Number} && error("Failed to transform $name to hardware, $name must be a vector Waveforms.\n If $name contains constant values use `BloqadeWaveforms.constant`.")
    if field isa Vector 
        foreach(field) do ele
            !(ele isa Waveform{F,T} where {F,T<:Real})  && error("Failed to transform $name to hardware, $name must be a Waveform.")
        end

        duration = field[1].duration
        foreach(field) do ele
            duration != ele.duration && error("Failed to transform $name to hardware, all Waveforms must have the same duration.")
        end

        return

    end
    (isnothing(field) || field isa Number || BloqadeExpr.is_time_function(field)) && error("Failed to transform $name to hardware, $name must be a Waveform.\n If $name is constant or zero use `BloqadeWaveforms.constant`.")

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
function pin_waveform_edges(wf::Waveform,name,
    max_slope::Real,
    begin_value::Real,
    end_value::Real)

    duration = wf.duration



    t_begin = if !isapprox(wf(0.0), begin_value;atol=eps(),rtol=√eps())
        @debug "During hardware transform: $name(t) start value is not $start_value. adding ramp(s) to fix endpoints."
        ramp_up =  (sign(wf(0.0)-begin_value)*max_slope)
        lin_ramp_begin = Waveform(t -> ramp_up .* t .+ begin_value, duration)
        find_zero(wf-lin_ramp_begin,(0.0,duration))
    else
        0.0
    end

    t_end = if !isapprox(wf(duration), end_value;atol=eps(),rtol=√eps())
        @debug "During hardware transform: $name(t) end value is not $end_value. adding ramp(s) to fix endpoints."
        ramp_down = (sign(end_value-wf(duration))*max_slope)
        lin_ramp_end = Waveform(t -> ramp_down .* (t.-duration) .+ end_value, duration)
        find_zero(wf-lin_ramp_end,(0.0,duration))
    else
        duration
    end

    # not the best solution but will 
    if t_begin > t_end
        t_begin = duration/2
        t_end = duration/2

        start_wf = linear_ramp(;
            duration=t_begin,
            start_value=begin_value,
            stop_value=wf(t_begin)
        )

        end_wf = linear_ramp(;
            duration=t_end,
            start_value=wf(t_end),
            stop_value=end_value
        )

        return append(start_wf,end_wf)        
    elseif t_begin > 0 && t_end < duration
        mid_wf = Waveform(t->wf.f(t.+t_begin),duration - t_begin - t_end)

        start_wf = linear_ramp(;
            duration=t_begin,
            start_value=begin_value,
            stop_value=wf(t_begin)
        )

        end_wf = linear_ramp(;
            duration=t_end,
            start_value=wf(t_end),
            stop_value=end_value
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

        end_wf = linear_ramp(;duration=t_end,
            start_value=wf(t_end),
            stop_value=end_value
        )

        return append(start_wf,end_wf)
    else
        return wf      
    end

end

# TODO: move this to BloqadeWaveforms
# does the SVD method to decompose local drives into an outer produce of masks and functions.
# function find_local_masks(values::Array{T,2};name::Symbol=:Waveform,ntrunc::Int=1, assert_truncation::Bool=false) where T
    
#     (l,w) = size(values)
    
#     # project out uniform amplitude waveforms
#     ones_norm = ones(w,1)./sqrt(w)
#     avg = (values*ones_norm)*transpose(ones_norm)
#     results = []
#     values = values - avg
#     avg = avg[:,1]

#     if ntrunc > 0
#         # use SVD to decompose the non-uniform waveforms
#         # if there are more than one singular value larger than 0 the remaining pattern
#         # can't be described as a local detuning mask times a single time-dependent function.
#         # here we truncate the waveform only keeping the first singular value.
#         u,s,v_T = svd(values)
#         contained_weight = sum(s[1:ntrunc])
        
#         remaining_weight = sum(s[(ntrunc+1):end])

#         if remaining_weight > eps()*length(s)*contained_weight
#             error_or_warn(assert_truncation,"Cannot decompose $name into product of masks and scalar functions.")
#         end

#         for i in 1:ntrunc

#             mask = if all(u[:,i] .<= 0) # make u all positive
#                 u[:,i] .*= -1
#                 (-s[i] .* v_T[:,i])

#             else
#                 s[i] .* v_T[:,i]
#             end

#             min = minimum(mask)
#             max = maximum(mask)

#             mask = (mask .- min) ./ (max - min)

#             f = u[:,i] .* (max - min)
#             avg .+= u[:,i] .* min

#             push!(results,(f,mask))
#         end
#     end

#     push!(results,(avg,ones(w)))
#     return results

# end
function clip_waveform(wf::Waveform{BloqadeWaveforms.PiecewiseLinear{T,I},T},name,min_value::T,max_value::T) where {T<:Real,I}
    @assert min_value < max_value

    for value in wf.f.values
        if value > max_value || value < min_value
            @debug "During hardware transform: $name(t) falls outside of hardware bounds, clipping to maximum/minimum."
            break
        end
    end

    return piecewise_linear(;clocks=wf.f.clocks,values=map(wf.f.values) do value
            return min(max_value,max(min_value,value))
        end
    )
end

function hardware_transform_Ω(Ω,device_capabilities::DeviceCapabilities)
    time_res = device_capabilities.rydberg.global_value.timeResolution
    min_step = device_capabilities.rydberg.global_value.timeDeltaMin
    rabi_res = device_capabilities.rydberg.global_value.rabiFrequencyResolution
    max_slope = device_capabilities.rydberg.global_value.rabiFrequencySlewRateMax
    max_value = device_capabilities.rydberg.global_value.rabiFrequencyMax
    min_value = device_capabilities.rydberg.global_value.rabiFrequencyMin

    check_waveform(Ω,:Ω)
    check_global(Ω,:Ω)
    warn_duration(time_res,Ω,:Ω)

    Ω = pin_waveform_edges(Ω,:Ω,max_slope,0,0)

    Ωt = if Ω isa PiecewiseLinearWaveform
        piecewise_linear(;
            clocks=set_resolution.(Ω.f.clocks,time_res),
            values=set_resolution.(Ω.f.values,rabi_res)
        )

    elseif Ω isa Waveform{F,T} where {F,T<:Real}
        
        Ω_interp = piecewise_linear_interpolate(Ω,max_slope=max_slope,min_step=min_step,atol=0)
        piecewise_linear(;
            clocks=set_resolution.(Ω_interp.f.clocks,time_res),
            values=set_resolution.(Ω_interp.f.values,rabi_res)
        )
    end

    Ωt = clip_waveform(Ωt,:Ω,min_value,max_value)
    return Ωt,norm_diff_durations(Ω,Ωt)
end

function hardware_transform_ϕ(ϕ,device_capabilities::DeviceCapabilities)
    time_res = device_capabilities.rydberg.global_value.timeResolution
    min_step = device_capabilities.rydberg.global_value.timeDeltaMin
    phase_res = device_capabilities.rydberg.global_value.phaseResolution
    max_slope = device_capabilities.rydberg.global_value.phaseSlewRateMax
    max_value = device_capabilities.rydberg.global_value.phaseMax
    min_value = device_capabilities.rydberg.global_value.phaseMin

    check_global(ϕ,:ϕ)
    check_waveform(ϕ,:ϕ)
    warn_duration(time_res,ϕ,:ϕ)

    ϕ = pin_waveform_edges(ϕ,:ϕ,max_slope,0.0,ϕ(ϕ.duration))

    ϕt = if ϕ isa PiecewiseLinearWaveform
        piecewise_linear(;
            clocks=set_resolution.(ϕ.f.clocks,time_res),
            values=set_resolution.(ϕ.f.values,phase_res)
        )        
    elseif ϕ isa Waveform{F,T} where {F,T<:Real}
        # arbitrary waveform must transform
        
        ϕ_interp = piecewise_linear_interpolate(ϕ,max_slope=max_slope,min_step=min_step,atol=0)
        piecewise_linear(;
            clocks=set_resolution.(ϕ_interp.f.clocks,time_res),
            values=set_resolution.(ϕ_interp.f.values,phase_res)
        )
    end

    ϕt = clip_waveform(ϕt,:ϕ,min_value,max_value)
    return ϕt,norm_diff_durations(ϕ,ϕt)
end

function hardware_transform_Δ(Δ,device_capabilities::DeviceCapabilities)

    time_res = device_capabilities.rydberg.global_value.timeResolution
    min_step = device_capabilities.rydberg.global_value.timeDeltaMin
    
    detune_res = device_capabilities.rydberg.global_value.detuningResolution
    local_res = device_capabilities.rydberg.local_value.localDetuningResolution
    common_detune_res = device_capabilities.rydberg.local_value.commonDetuningResolution
    max_slope = device_capabilities.rydberg.global_value.detuningSlewRateMax
    max_value = device_capabilities.rydberg.global_value.detuningMax
    min_value = device_capabilities.rydberg.global_value.detuningMin


    check_waveform(Δ,:Δ)
    check_global(Δ,:Δ)
    warn_duration(time_res,Δ,:Δ)
    


    if Δ isa PiecewiseLinearWaveform
        Δt = piecewise_linear(;
            clocks=set_resolution.(Δ.f.clocks,time_res),
            values=set_resolution.(Δ.f.values,detune_res)
        )
        Δ_mask = (Δ=Δt,δ=nothing,Δi=1.0)
        Δt = clip_waveform(Δt,:Δ,min_value,max_value)
        return Δt,norm_diff_durations(Δ,Δt),Δ_mask

    elseif Δ isa Waveform{F,T} where {F,T<:Real}
        # arbitrary waveform must transform
        
        Δ_interp = piecewise_linear_interpolate(Δ,max_slope=max_slope,min_step=min_step,atol=0)
        Δt = piecewise_linear(;
            clocks=set_resolution.(Δ_interp.f.clocks,time_res),
            values=set_resolution.(Δ_interp.f.values,detune_res)
        )

        Δt = clip_waveform(Δt,:Δ,min_value,max_value)
        Δ_mask = (Δ=Δt,δ=nothing,Δi=1.0)

        return Δt,norm_diff_durations(Δ,Δt),Δ_mask
    elseif Δ isa Vector
        throw(NotImplementedError("Not implemented."))

        """
        duration = Δ[1].duration
        nsteps = Int((duration+duration%min_step)÷min_step)
        clocks = collect(LinRange(0.0,duration,nsteps))
        clocks = set_resolution.(clocks,time_res)

        values = zeros(length(clocks),length(Δ))

        @inbounds for (j,δ) in enumerate(Δ)
            values[:,j] = δ.(clocks)
        end

        ((δ_values,Δi),(Δ_values,_)) = find_local_masks(values;name=:Δ)

        Δ_values = set_resolution.(Δ_values,detune_res)
        δ_values = set_resolution.(δ_values,common_detune_res)
        Δi = set_resolution.(Δi,local_res)

        Δ_mask = (
            Δ=piecewise_linear(;clocks=clocks,values=Δ_values),
            δ=piecewise_linear(;clocks=clocks,values=Δ_values),
            Δi=Δi
        )
        Δt = [piecewise_linear(;clocks=clocks,values=Δ_values+δi.*δ_values) for δi in Δi]
        error = [norm_diff_durations(δ,δt) for (δ,δt) in zip(Δ,Δt)]

        return Δt,error,Δ_mask
        """
    end


end

function hardware_transform_parse(H::BloqadeExpr.RydbergHamiltonian,device_capabilities::DeviceCapabilities)
    (atoms,ϕ,Ω,Δ) = get_rydberg_params(H)

    ϕ,ϕ_error = hardware_transform_ϕ(ϕ,device_capabilities)
    Ω,Ω_error = hardware_transform_Ω(Ω,device_capabilities)
    Δ,Δ_error,Δ_mask = hardware_transform_Δ(Δ,device_capabilities)

    pos_resolution = device_capabilities.lattice.geometry.positionResolution
    new_atoms = [set_resolution.(pos,pos_resolution) for pos in atoms]

    mse_atoms = sum(√sum((a .- b) .^ 2) for (a,b) in zip(new_atoms,atoms))/length(atoms)

    info = (ϕ=ϕ_error,Ω=Ω_error,Δ=Δ_error,Δ_mask=Δ_mask,mse_atoms=mse_atoms)

    return (new_atoms,ϕ,Ω,Δ,info)
end

# public API exposed here: 
function hardware_transform(H::BloqadeExpr.RydbergHamiltonian;device_capabilities::DeviceCapabilities=get_device_capabilities())
    atoms,ϕ,Ω,Δ,info = hardware_transform_parse(H,device_capabilities)
    hardware_H = rydberg_h(atoms,ϕ=ϕ,Ω=Ω,Δ=Δ)

    return hardware_H,info
end
