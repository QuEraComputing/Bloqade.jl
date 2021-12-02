struct CompositeWaveform{N, T, Waveforms <: Tuple} <: AbstractWaveform
    waveforms::Waveforms
    durations::NTuple{N, T}
    checkpoints::NTuple{N, T}
end

function CompositeWaveform(waveforms::AbstractWaveform...)
    durations = ntuple(duration, waveforms)
    return CompositeWaveform(waveforms, durations, accumulate(+, durations))
end

duration(waveform::CompositeWaveform) = waveform.checkpoints[end]

# NOTE: change this to a generated function if the performance matters
# in the future
function (waveform::CompositeWaveform)(t::Real, offset::Real=zero(t))
    clock = t - offset
    local idx
    for idx in 1:length(waveform.checkpoints)
        if clock ≤ waveform.checkpoints[idx]
            break
        end
    end
    return waveform.waveforms[idx](t, offset + waveform.checkpoints[idx-1])
end
