"""
abstract type AbstractWaveform end

Waveforms are callable structs that has a [`duration`](@ref)
and a callable method `(::WaveForm)(t::Real, offset::Real=zero(t))`.

Where the argument `t` is the clock and `offset` is an offset to
the clock for composition.
"""
abstract type AbstractWaveform end

function duration end

function assert_clock(t::Real, duration::Real, offset::Real)
    min, max = offset, duration + offset
    msg = "clock out of range $min ≤ t ≤ $max, got clock $t"
    min ≤ t ≤ max || throw(BoundsError(msg))
    return
end
