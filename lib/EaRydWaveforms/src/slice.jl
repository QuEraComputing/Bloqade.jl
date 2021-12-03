struct SlicedWaveform{T, Waveform, IntervalType <: Interval{<:Real, Closed, Closed}} <: AbstractWaveform
    waveform::Waveform
    interval::IntervalType
end

function SlicedWaveform(waveform, interval::Interval{<:Real, Closed, Closed})
    return SlicedWaveform{eltype(waveform), typeof(waveform), typeof(interval)}(waveform, interval)
end

Base.eltype(wf::SlicedWaveform) = eltype(wf.waveform)
duration(wf::SlicedWaveform) = span(wf.interval)

function (wf::SlicedWaveform)(t::Real, offset::Real=zero(t))
    assert_clock(t, duration(wf), offset)
    return wf.waveform(t, offset - wf.interval.first)
end

function Base.getindex(wf::AbstractWaveform, slice::Interval)
    0 ≤ slice.first ≤ duration(wf) && 0 ≤ slice.last ≤ duration(wf) ||
        throw(ArgumentError(
            "slice out of range, expect 0 ≤ slice ≤ $(duration(wf)) " *
            "got $slice"))
    return SlicedWaveform(wf, slice)
end
