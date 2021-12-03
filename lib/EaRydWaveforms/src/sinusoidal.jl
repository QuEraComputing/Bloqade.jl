Base.@kwdef struct SinusoidalWaveform{T <: Real} <: AbstractWaveform
    duration::T
    amplitude::T = 1.0
end

Base.eltype(::SinusoidalWaveform{T}) where T = T
duration(waveform::SinusoidalWaveform) = waveform.duration

function (waveform::SinusoidalWaveform{T})(t::Real, offset::Real=zero(t)) where T
    assert_clock(t, duration(waveform), offset)
    return waveform.amplitude * sin(T(t - offset))
end
