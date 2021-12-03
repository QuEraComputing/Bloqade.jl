struct FunctionWaveform{F, T <: Real} <: AbstractWaveform
    f::F
    duration::T
end

function FunctionWaveform(f; duration::Real)
    return FunctionWaveform(f, duration)
end

duration(wf::FunctionWaveform) = wf.duration
Base.eltype(::FunctionWaveform{<:Any, T}) where T = T

function (wf::FunctionWaveform)(t::Real, offset::Real=zero(t))
    assert_clock(t, duration(wf), offset)
    return wf.f(t - offset)
end
