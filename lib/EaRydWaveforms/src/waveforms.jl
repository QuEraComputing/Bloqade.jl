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

function Base.show(io::IO, mime::MIME"text/plain", waveform::AbstractWaveform)
    name = nameof(typeof(waveform))
    xs = range(0, duration(waveform), 100)
    plt = lineplot(
        xs, waveform.(xs);
        title=string(name, "{", eltype(waveform), "}"),
        # TODO: decide the unit?
        xlabel="clock (μs)",
        ylabel="value (rad/µs)",
        compact=true,
    )
    return show(io, mime, plt)
end
