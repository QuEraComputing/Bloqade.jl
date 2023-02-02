# this function is required to calculate difference
# between waveforms which have different durations.
# the behavior is to pad the shorter waveform at the
# end with 0. This is useful when the user passes
# in a duration that doesn't neatly fit into the
# time resolution of the device. In which case 
# this is required to check the error.
function norm_diff_durations(A::Waveform,B::Waveform)
    A.duration==B.duration && return norm(A-B)

    if A.duration > B.duration
        T_diff = A.duration - B.duration
        B = append(B,constant(;duration=T_diff,value=0))
        return norm(A-B)
    else
        T_diff = B.duration - A.duration
        A = append(A,constant(;duration=T_diff,value=0))
        return norm(A-B)
    end
end

# simple function, round number to the nearest multiple of `res`
function set_resolution(x::Real,res::Real)
    @assert res > 0
    return round(Int(round(x / res)) * res,sigdigits=14)
end

# throws error for non-scalar fields
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
        @debug "waveform $name duration will be rounded during hardware transformation."
    end
end

# TODO: move this to BloqadeWaveforms
# function inserts ramps at the beginning and end to set correct initial and final values
# two methods, one for general waveform, other is for specifically PWL waveforms
function pin_waveform_edges(wf::Waveform,name,
    max_slope::Real,
    begin_value::Real,
    end_value::Real)

    
    duration = wf.duration

    wf_begin = wf(zero(duration))
    wf_end = wf(duration)

    
    # 
    t_begin = if !isapprox(wf_begin, begin_value;atol=eps(),rtol=√eps())
        @debug "During hardware transform: $name(t) initial value is not $begin_value. adding ramp to fix endpoints."
        ramp_up =  (sign(wf_begin-begin_value)*max_slope)
        lin_ramp_begin = Waveform(t -> ramp_up .* t .+ begin_value, duration)
        find_zero(wf-lin_ramp_begin,(zero(duration),duration))
    else
        zero(duration)
    end

    t_end = if !isapprox(wf_end, end_value;atol=eps(),rtol=√eps())
        @debug "During hardware transform: $name(t) end value is not $end_value. adding ramp to fix endpoints."
        ramp_down = (sign(end_value-wf_end)*max_slope)
        lin_ramp_end = Waveform(t -> ramp_down * (t-duration) + end_value, duration)
        find_zero(wf-lin_ramp_end,(zero(duration),duration))
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
        mid_wf = wf[t_begin..t_end]

        start_wf = linear_ramp(;
            duration=t_begin,
            start_value=begin_value,
            stop_value=wf(t_begin)
        )

        end_wf = linear_ramp(;
            duration=duration-t_end,
            start_value=wf(t_end),
            stop_value=end_value
        )


        return append(start_wf,mid_wf,end_wf)
    elseif t_begin > 0 
        end_wf = wf[t_begin..duration]
        start_wf = linear_ramp(;duration=t_begin,
            start_value=begin_value,
            stop_value=wf(t_begin)
        )

        return append(start_wf,end_wf)
    elseif t_end < duration
        start_wf = wf[0..t_end]

        end_wf = linear_ramp(;duration=duration-t_end,
            start_value=wf(t_end),
            stop_value=end_value
        )

        return append(start_wf,end_wf)
    else
        return wf      
    end

end

#=
function find_local_masks(values::Array{T,2};name::Symbol=:Waveform,ntrunc::Int=1, assert_truncation::Bool=false) where T
    
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
        contained_weight = sum(s[1:ntrunc])
        
        remaining_weight = sum(s[(ntrunc+1):end])

        if remaining_weight > eps()*length(s)*contained_weight
            error_or_warn(assert_truncation,"Cannot decompose $name into product of masks and scalar functions.")
        end

        for i in 1:ntrunc

            mask = if all(u[:,i] .<= 0) # make u all positive
                u[:,i] .*= -1
                (-s[i] .* v_T[:,i])

            else
                s[i] .* v_T[:,i]
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
=#

function clip_waveform(wf::Waveform{BloqadeWaveforms.PiecewiseLinear{T,I},T},name,min_value::T,max_value::T) where {T<:Real,I}
    #=
    FROM JULIA DOCS:
    !!! warning An assert might be disabled at various optimization levels. 
    Assert should therefore only be used as a debugging tool and not used for 
    authentication verification (e.g., verifying passwords), 
    nor should side effects needed for the function to work correctly be used inside of asserts.
    =#
    # ensure that the maximum value of waveform is greater than minimum value
    @assert min_value < max_value

    # for each value, if a value is greater than the maximum permitted or smaller than the smallest permitted, 
    # issue notice to user that clipping must occur
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

function clip_waveform(wf::Waveform{BloqadeWaveforms.PiecewiseConstant{T},T},name,min_value::T,max_value::T) where {T<:Real}
    @assert min_value < max_value

    for value in wf.f.values
        if value > max_value || value < min_value
            @debug "During hardware transform: $name(t) falls outside of hardware bounds, clipping to maximum/minimum."
            break
        end
    end

    return piecewise_constant(;clocks=wf.f.clocks,values=map(wf.f.values) do value
            return min(max_value,max(min_value,value))
        end
    )
end


function hardware_transform_parse(h::BloqadeExpr.RydbergHamiltonian,device_capabilities::DeviceCapabilities)
    (atoms,ϕ,Ω,Δ) = get_rydberg_params(h)
    Δ_mask = (Δ=Δ,δ=nothing,Δi=1.0)

    ϕ,ϕ_error = hardware_transform_ϕ(ϕ,device_capabilities)
    Ω,Ω_error = hardware_transform_Ω(Ω,device_capabilities)
    Δ,Δ_error = hardware_transform_Δ(Δ,device_capabilities)
    atoms,mse_atoms = hardware_transform_atoms(atoms,device_capabilities)

    info = HardwareTransformInfo(ϕ_error=ϕ_error,Ω_error=Ω_error,Δ_error=Δ_error,Δ_mask=Δ_mask,mse_atoms=mse_atoms)

    return (atoms,ϕ,Ω,Δ,info)
end

# public API exposed here: 
"""
    hardware_transform_Ω(Ω,device_capabilities::DeviceCapabilities=get_device_capabilities())

Given the `device_capabilities` of the machine and a Rabi frequency (Ω) [`Waveform`](@ref),
return a transformed Ω capable of being implemented by the machine along with
the error between the original (``A``) and transformed (``B``)
waveforms calculated as ``\\Vert A - B\\Vert_1``. If the waveform durations are different,
the shorter waveform is padded with zeros for values to make the durations equal in error calculation.

# Logs/Warnings/Exceptions

Exceptions are thrown if Ω is:
- Not of type [`BloqadWaveforms.Waveform`](@ref)
- Not a global drive (e.g.: Vector of Waveforms, localized Ω is not supported)
- the maximum slope allowed for the waveform from `device_capabilities` is set to infinity 
- the minimum time step allowed for the waveform from `device_capabilities` is set to zero


Debug logs are issued if the following are encountered in Ω:
- duration may be rounded due to time resolution from `device_capabilities`
- the initial waveform does not start/end in zero for its value
- the values in the waveform exceed `device_capabilities` supported values, and must be clipped 

# Examples
```jldoctest; setup:=(using BloqadeWaveforms)
julia> wf = sinusoidal(duration=2, amplitude=1.3*π);

julia> hardware_transform_Ω(wf)
(Waveform(_, 2.0), 2.632451578170084)
```
"""
function hardware_transform_Ω(Ω,device_capabilities::DeviceCapabilities=get_device_capabilities())
    time_res = device_capabilities.rydberg.global_value.time_resolution
    min_step = device_capabilities.rydberg.global_value.time_delta_min
    rabi_res = device_capabilities.rydberg.global_value.rabi_frequency_resolution
    max_slope = device_capabilities.rydberg.global_value.rabi_frequency_slew_rate_max
    max_value = device_capabilities.rydberg.global_value.rabi_frequency_max
    min_value = device_capabilities.rydberg.global_value.rabi_frequency_min

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

"""
    hardware_transform_ϕ(ϕ,device_capabilities::DeviceCapabilities=get_device_capabilities())

Given the `device_capabilities` of the machine and a phase [`Waveform`](@ref),
return a transformed ϕ capable of being implemented by the machine along with
the error between the original (``A``) and transformed (``B``)
waveforms calculated as ``\\Vert A - B\\Vert_1``. If the waveform durations are different,
the shorter waveform is padded with zeros for values to make the durations equal in error calculation.

# Logs/Warnings/Exceptions

Exceptions are thrown if ϕ is:
- Not of type [`BloqadWaveforms.Waveform`](@ref)
- Not a global drive (e.g.: Vector of Waveforms, localized ϕ is not supported)
- the maximum slope allowed for the waveform from `device_capabilities` is set to infinity 
- the minimum time step allowed for the waveform from `device_capabilities` is set to zero


Debug logs are issued if the following are encountered in ϕ:
- duration may be rounded due to time resolution from `device_capabilities`
- the initial waveform does not start/end in zero for its value
- the values in the waveform exceed `device_capabilities` supported values, and must be clipped 

# Examples
```jldoctest; setup:=(using BloqadeWaveforms)
julia> wf = sinusoidal(duration=2, amplitude=1.3*π);

julia> hardware_transform_ϕ(wf)
(Waveform(_, 2.0), 0.5386117854062276)
```
"""
function hardware_transform_ϕ(ϕ,device_capabilities::DeviceCapabilities=get_device_capabilities())

    time_res = device_capabilities.rydberg.global_value.time_resolution
    min_step = device_capabilities.rydberg.global_value.time_delta_min
    phase_res = device_capabilities.rydberg.global_value.phase_resolution
    max_value = device_capabilities.rydberg.global_value.phase_max
    min_value = device_capabilities.rydberg.global_value.phase_min

    # Can't feed in array of waveforms for ϕ
    check_global(ϕ,:ϕ)
    # check that it actually is a waveform (ex: an int/float constant can't be a waveform)
    check_waveform(ϕ,:ϕ)
    # warn if duration will be rounded to time resolution
    warn_duration(time_res,ϕ,:ϕ)

    ϕt = if ϕ isa PiecewiseConstantWaveform
        piecewise_constant(;
            clocks=set_resolution.(ϕ.f.clocks,time_res),
            values=set_resolution.(ϕ.f.values,phase_res)
        )        

    elseif ϕ isa Waveform{F,T} where {F,T<:Real} # handle the more general case
        # arbitrary waveform must transform
        ϕ_interp = piecewise_constant_interpolate(ϕ,min_step=min_step,atol=0)
        piecewise_constant(;
            clocks=set_resolution.(ϕ_interp.f.clocks,time_res),
            values=set_resolution.(ϕ_interp.f.values,phase_res)
        )
    end 
    # simply clip value to make waveform consistent
    if !isapprox(ϕt(zero(ϕt.duration)),zero(ϕt.duration);atol=eps()) 
        ϕt.f.values[1] = zero(ϕt.duration)
    end

    ϕt = clip_waveform(ϕt,:ϕ,min_value,max_value)
    return ϕt,norm_diff_durations(ϕ,ϕt)
end

"""
    hardware_transform_Δ(Δ,device_capabilities::DeviceCapabilities=get_device_capabilities())

Given the `device_capabilities` of the machine and a detuning waveform Δ,
return a transformed Δ capable of being implemented by the machine along with
the error between the original (``A``) and transformed (``B``)
waveforms calculated as ``\\Vert A - B\\Vert_1``. If the waveform durations are different,
the shorter waveform is padded with zeros for values to make the durations equal in error calculation.

# Logs/Warnings/Exceptions

Exceptions are thrown if Δ is:
- Not of type [`BloqadWaveforms.Waveform`](@ref)
- Not a global drive (e.g. Vector of Waveforms, localized Δ is not supported)
- the maximum slope allowed for the waveform from `device_capabilities` is set to infinity 
- the minimum time step allowed for the waveform from `device_capabilities` is set to zero


Debug logs are issued if the following are encountered in Δ:
- duration may be rounded due to time resolution from `device_capabilities`
- the initial waveform does not start/end in zero for its value
- the values in the waveform exceed `device_capabilities` supported values, and must be clipped 

# Examples
```jldoctest; setup:=(using BloqadeWaveforms)
julia> wf = sinusoidal(duration=2, amplitude=1.3*π);

julia> hardware_transform_Δ(wf)
(Waveform(_, 2.0), 0.06492452289703464)
```
"""
function hardware_transform_Δ(Δ,device_capabilities::DeviceCapabilities=get_device_capabilities())

    time_res = device_capabilities.rydberg.global_value.time_resolution
    min_step = device_capabilities.rydberg.global_value.time_delta_min
    
    detune_res = device_capabilities.rydberg.global_value.detuning_resolution
    # local_res = device_capabilities.rydberg.local_value.local_detuning_resolution
    # common_detune_res = device_capabilities.rydberg.local_value.common_detuning_resolution
    max_slope = device_capabilities.rydberg.global_value.detuning_slew_rate_max
    max_value = device_capabilities.rydberg.global_value.detuning_max
    min_value = device_capabilities.rydberg.global_value.detuning_min


    check_waveform(Δ,:Δ)
    check_global(Δ,:Δ)
    warn_duration(time_res,Δ,:Δ)
    


    if Δ isa PiecewiseLinearWaveform
        Δt = piecewise_linear(;
            clocks=set_resolution.(Δ.f.clocks,time_res),
            values=set_resolution.(Δ.f.values,detune_res)
        )
        Δt = clip_waveform(Δt,:Δ,min_value,max_value)
        return Δt,norm_diff_durations(Δ,Δt)

    elseif Δ isa Waveform{F,T} where {F,T<:Real}
        # arbitrary waveform must transform
        
        Δ_interp = piecewise_linear_interpolate(Δ,max_slope=max_slope,min_step=min_step,atol=0)
        Δt = piecewise_linear(;
            clocks=set_resolution.(Δ_interp.f.clocks,time_res),
            values=set_resolution.(Δ_interp.f.values,detune_res)
        )

        Δt = clip_waveform(Δt,:Δ,min_value,max_value)

        return Δt,norm_diff_durations(Δ,Δt)
    elseif Δ isa Vector
        error("Local detuning not implemented in Schema.")
        #=
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
        =#
    end


end

"""
    hardware_transform_atoms(atoms,device_capabilities::DeviceCapabilities=get_device_capabilities())

Given the constraints of the hardware from `device_capabilities` (specifically the position resolution)
and an iterable containing the atom positions `atoms`, returns a Tuple containing the adjusted 
atom positions the machine is capable of resolving and the mean squared error between the 
desired atom positions and newly generated ones.

Note that other constraints such as the maximum width, height, and minimum supported spacings
are not taken into account in adjusting the atoms. This may result in the [`validation`](@ref) 
function failing and requiring user intervention to modify the atom positions such that 
they satisfy the other constraints.

# Examples 
```jldoctest;
julia> atom_positions = ((1.12,), (2.01,), (3.01,));

julia> hardware_transform_atoms(atom_positions) # by default, calls get_device_capabilities()
([(1.1,), (2.0,), (3.0,)], 0.013333333333333197)
```
"""
function hardware_transform_atoms(atoms,device_capabilities::DeviceCapabilities=get_device_capabilities())

    pos_resolution = device_capabilities.lattice.geometry.position_resolution

    new_atoms = [set_resolution.(pos,pos_resolution) for pos in atoms]

    mse_atoms = sum(√sum((x.-y).^2) for (x,y) in zip(new_atoms,atoms))/length(new_atoms)

    return new_atoms,mse_atoms

end

"""
    hardware_transform(h::BloqadeExpr.RydbergHamiltonian;device_capabilities::DeviceCapabilities=get_device_capabilities())

Given the constraints of the hardware via `device_capabilities`, transforms `h` into one the
machine is capable of executing as well as:
* The mean squared error between original positions of the atoms and the transformed ones
* The 1-norm of the difference between the original and transformed waveforms
which are all stored in a [`HardwareTransformInfo`](@ref) struct.

Note that not all atom position constraints are accounted for, such as the maximum lattice width, lattice height, 
and minimum supported spacings. Only position resolution is automatically accounted for.
This may result in the [`validation`](@ref) function failing and requiring user intervention to modify the atom 
positions such that they satisfy the other constraints.

# Logs/Warnings/Exceptions

Debug logs are *always* emitted containing the error (defined as the 1-norm of the difference between
the original waveform and the transformed waveform) across all waveforms (Ω, Δ, ϕ) as well as the 
Mean Squared Error between the original atom positions and the adjusted atom positions.

The debug logs/warnings from constituent functions [`hardware_transform_Ω`](@ref),
[`hardware_transform_Δ`](@ref), [`hardware_transform_ϕ`](@ref) are also emitted should the 
waveforms in `h` cause them to.

See also [`hardware_transform_atoms`](@ref), [`hardware_transform_Ω`](@ref),
[`hardware_transform_Δ`](@ref), [`hardware_transform_ϕ`](@ref).

# Examples
```jldoctest; setup:=(using BloqadeExpr, BloqadeLattices)
julia> atom_positions = AtomList([(1.12,), (2.01,), (3.01,)]);

julia> Δ = Ω = ϕ = sinusoidal(duration=2, amplitude=1.3*π);

julia> h = rydberg_h(atom_positions; Ω=Ω,Δ=Δ,ϕ=ϕ)
nqubits: 3
+
├─ [+] ∑ 2π ⋅ 8.627e5.0/|x_i-x_j|^6 n_i n_j
├─ [+] Ω(t) ⋅∑ e^{ϕ(t) ⋅ im} |0⟩⟨1| + e^{-ϕ(t) ⋅ im} |1⟩⟨0|
└─ [-] Δ(t) ⋅ ∑ n_i


julia> hardware_transform(h)
(nqubits: 3
+
├─ [+] ∑ 2π ⋅ 8.627e5.0/|x_i-x_j|^6 n_i n_j
├─ [+] Ω(t) ⋅∑ e^{ϕ(t) ⋅ im} |0⟩⟨1| + e^{-ϕ(t) ⋅ im} |1⟩⟨0|
└─ [-] Δ(t) ⋅ ∑ n_i
, BloqadeSchema.HardwareTransformInfo(0.5386117854062276, 2.632451578170084, 0.06492452289703464, (Δ = Waveform(_, 2), δ = nothing, Δi = 1.0), 0.013333333333333197))
```
"""
function hardware_transform(h::BloqadeExpr.RydbergHamiltonian;device_capabilities::DeviceCapabilities=get_device_capabilities())
    # gets transformed versions of waveforms, then creates new hamiltonian    
    atoms,ϕ,Ω,Δ,info = hardware_transform_parse(h,device_capabilities)
    hardware_h = rydberg_h(atoms,ϕ=ϕ,Ω=Ω,Δ=Δ)

    @debug "Hardware transform report: after linear interpolation ∫dt |ϕ(t)-ϕ_hw(t)| = $(info.ϕ_error) rad⋅μs"
    @debug "Hardware transform report: after linear interpolation ∫dt |Ω(t)-Ω_hw(t)| = $(info.Ω_error) rad"
    @debug "Hardware transform report: after linear interpolation ∫dt |Δ(t)-Δ_hw(t)| = $(info.Δ_error) rad"
    @debug "Hardware transform report: mean deviation after rounding positions $(info.mse_atoms) μm"

    return hardware_h,info
end
