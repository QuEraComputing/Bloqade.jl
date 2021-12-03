"""
abstract type AbstractWaveform end

Abstract type for waveforms.

# Interfaces

Waveforms support the following interface:

- [`duration`](@ref): the duration of waveform.
- a callable method `(::WaveForm)(t::Real, offset::Real=zero(t))` returns
    the value of the waveform at time `t` with an given `offset` of the time.
- `Base.eltype`: returns the type of each value.
- [`sample_values`](@ref): sample the values of waveform.
"""
abstract type AbstractWaveform end

function duration end

function assert_clock(t::Real, duration::Real, offset::Real)
    min, max = offset, duration + offset
    msg = "clock out of range $min ≤ t ≤ $max, got clock $t"
    min ≤ t ≤ max || throw(ArgumentError(msg))
    return
end

function sample_values(waveform::AbstractWaveform, dt::Real=1e-3)
    clocks = range(zero(dt), duration(waveform); step=dt)
    return map(waveform, clocks)
end

function Base.show(io::IO, waveform::AbstractWaveform)
    summary(io, waveform)
    print(io, "(...)")
end

function Base.show(io::IO, mime::MIME"text/plain", waveform::AbstractWaveform)
    name = nameof(typeof(waveform))
    xs = 0:1e-3:duration(waveform)
    plt = lineplot(
        xs, sample_values(waveform);
        title=summary(waveform),
        # TODO: decide the unit?
        xlabel="clock (μs)",
        ylabel="value (rad/µs)",
        compact=true,
    )
    return show(io, mime, plt)
end
