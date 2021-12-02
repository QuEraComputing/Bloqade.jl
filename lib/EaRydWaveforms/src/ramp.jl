Base.@kwdef struct RampWaveform{T <: Real} <: AbstractWaveform
    duration::T
    start::T # start value
    stop::T # stop value
end

duration(waveform::RampWaveform) = waveform.duration

function (waveform::RampWaveform)(t::Real, offset::Real=zero(t))
    assert_clock(t, duration(waveform), offset)
    return (waveform.stop - waveform.start) / duration(waveform) * (t - offset)
end
