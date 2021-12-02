struct CompositeWaveform{N, T, Waveforms <: Tuple} <: AbstractWaveform
    waveforms::Waveforms
    durations::NTuple{N, T}
    checkpoints::NTuple{N, T}
end

function CompositeWaveform(waveforms::AbstractWaveform...)
    durations = map(duration, waveforms)
    return CompositeWaveform(waveforms, durations, accumulate(+, durations))
end

Base.eltype(waveform::CompositeWaveform) = eltype(waveform.waveforms[1])
Base.getindex(waveform::CompositeWaveform, idx::Int) = waveform.waveforms[idx]
Base.length(waveform::CompositeWaveform) = length(waveform.waveforms)

duration(waveform::CompositeWaveform) = waveform.checkpoints[end]

function sample_values(waveform::CompositeWaveform, dt::Real)
    values = typeof(dt)[]
    for wf in waveform.waveforms
        for clock in 0:dt:duration(wf)-dt
            push!(values, wf(clock))
        end
    end

    # push last point to match the inteval
    # length
    wf = last(waveform.waveforms)
    push!(values, wf(duration(wf)))
    return values
end

# NOTE: change this to a generated function if the performance matters
# in the future
function (waveform::CompositeWaveform)(t::Real, offset::Real=zero(t))
    clock = t - offset
    local idx
    for waveform_idx in 1:length(waveform.checkpoints)
        if clock ≤ waveform.checkpoints[waveform_idx]
            idx = waveform_idx
            break
        end
    end

    offset = isone(idx) ? offset : offset + waveform.checkpoints[idx-1]
    return waveform.waveforms[idx](t, offset)
end
