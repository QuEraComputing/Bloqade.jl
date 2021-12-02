struct SinusoidalWaveform{T <: Real} <: AbstractWaveform
    duration::T
end

Base.eltype(::SinusoidalWaveform{T}) where T = T
duration(waveform::SinusoidalWaveform) = waveform.duration

function (waveform::SinusoidalWaveform)(t::Real, offset::Real=zero(t))
    assert_clock(t, duration(waveform), offset)
    return sin(t - offset)
end
