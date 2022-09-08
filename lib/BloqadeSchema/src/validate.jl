
get_device_capabilities(kw...) = DeviceCapabilities(kw...)

function validate_lattice(positions,lattice_capabilities::LatticeCapabilities) end

function get_field_capabilities(field::Symbol,rydberg::RydbergCapabilities)

    if field == :Ω
        return (
            min_time_step = rydberg.global.timeDeltaMin,
            time_resolution = rydberg.global.timeResolution,
            max_time = rydberg.global.timeMax,
            max_value = rydberg.global.rabiFrequencyMax,
            min_value = rydberg.global.rabiFrequencyMin,
            max_slope = rydberg.global.rabiFrequencySlewRateMax,
            value_resolution = rydberg.global.rabiFrequencyResolution
            )
    elseif field == :ϕ
        return (
            min_time_step = rydberg.global.timeDeltaMin,
            time_resolution = rydberg.global.timeResolution,
            max_time = rydberg.global.timeMax,
            max_value = rydberg.global.phaseMax,
            min_value = rydberg.global.phaseMin,
            # max_slope = rydberg.global.phaseSlewRateMax,
            value_resolution = rydberg.global.phaseResolution
            )
    elseif field == :Δ
        return (
            min_time_step = rydberg.global.timeDeltaMin,
            time_resolution = rydberg.global.timeResolution,
            max_time = rydberg.global.timeMax,
            max_value = rydberg.global.detuningMax,
            min_value = rydberg.global.detuningMin,
            max_slope = rydberg.global.detuningSlewRateMax,
            value_resolution = rydberg.global.detuningResolution
        )
    else
        error("Symbol $(field) not recognized. Must be :Ω, :ϕ, or :Δ.")
    end

end





function validate_Ω(wf::Waveform{PiecewiseLinear{T,Interp},T},field::Symbol,rydberg::RydbergCapabilities) where {T<:Real,Interp}
    capabilities = get_field_capabilities(:Ω,rydberg)


    # validate duration
    # validate min_step
    # validate max slope
    # validate min value
    # validate max value
    # validate begin value
    # validate end value
end
validate_Ω(::Nothing,::RydbergCapabilities) = nothing
validate_Ω(::Any,::RydbergCapabilities) = error("Schema only supports piecewise linear and constant waveforms for Rydberg fields.")


function validate_Δ(wf::Waveform{PiecewiseLinear{T,Interp},T},field::Symbol,rydberg::RydbergCapabilities) where {T<:Real,Interp}
    capabilities = get_field_capabilities(:Δ,rydberg)

    # validate duration
    # validate min_step
    # validate max slope
    # validate min value
    # validate max value

end
validate_Δ(::Nothing,::RydbergCapabilities) = nothing
validate_Δ(::Any,::RydbergCapabilities) = error("Schema only supports piecewise linear and constant waveforms for Rydberg fields.")


function validate_ϕ(wf::Waveform{PiecewiseConstant{T},T}) where {T<:Real} 
    capabilities = get_field_capabilities(:ϕ,rydberg)

    # validate duration
    # validate min_step
    # validate min value
    # validate max value
end
validate_ϕ(::Nothing,::RydbergCapabilities) = nothing
validate_ϕ(::Any,::RydbergCapabilities) = error("Schema only supports piecewise linear and constant waveforms for Rydberg fields.")


function mult_by_two(Ω)
    isnothing(Ω) && return
    if BloqadeExpr.is_const_param(Ω)
        return Ω .* 2
    end

    return if Ω isa Vector
        map(Ω) do Ω_i
            return Ω_i.f
        end
    else
        Ω.f
    end

end

function validate_waveform(::Nothing;kw...) end

function validate(H::BloqadeExpr.RydbergHamiltonian;kw...)
    device_capabilities = get_device_capabilities(kw...)

    ϕ = nothing
    Ω = nothing
    Δ = nothing


    validate_lattice(H.rydberg_term.atoms,device_capabilities.lattice)

    if H.rabi_term isa BloqadeExpr.SumOfX
        Ω = mult_by_two(H.rabi_term.Ω)  
    elseif H.rabi_term isa BloqadeExpr.SumOfXPhase
        Ω = mult_by_two(H.rabi_term.Ω)
        ϕ = H.rabi_term.ϕ

    end

    if H.detuning_term isa BloqadeExpr.SumOfN
        Δ = H.detuning_term.Δ
    end

    validate_Ω(Ω,device_capabilities.rydberg)
    validate_ϕ(ϕ,device_capabilities.rydberg)
    validate_Δ(Δ,device_capabilities.rydberg)

end
