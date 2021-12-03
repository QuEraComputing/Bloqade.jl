Base.@kwdef struct ConstantWaveform{T <: Real} <: AbstractWaveform
    amplitude::T
    duration::T
end

duration(wf::ConstantWaveform) = wf.duration
Base.eltype(::ConstantWaveform{T}) where T = T
function (wf::ConstantWaveform)(t::Real, offset::Real=zero(t))
    assert_clock(t, wf.duration, offset)
    return wf.amplitude
end
