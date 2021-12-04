using Intervals

struct Waveform{F, T <: Interval{<:Real, Closed, Closed}}
    f::F
    interval::T
end

Waveform(f; start::Real=0, stop::Real) = Waveform(f, start..stop)
Base.eltype(wf::Waveform) = eltype(wf.interval)

function Base.getindex(wf::Waveform, slice::Interval{<:Real, Closed, Closed})
    issubset(slice, wf.interval) || throw(ArgumentError("slice is not in $(wf.interval) range, got $slice"))
    return Waveform(wf.f, slice)
end

function (wf::Waveform)(t::Real)
    t in wf.interval || throw(ArgumentError("t is not in $(wf.interval) range, got $t"))
    return wf.f(t)
end

function sample_clock(wf::Waveform, dt::Real=1e-3)
    return wf.interval.first:dt:wf.interval.last
end

function sample_values(wf::Waveform, dt::Real=1e-3)
    return [wf(t) for t in sample_clock(wf, dt)]
end

Base.show(io::IO, wf::Waveform) = print(io, "Waveform(", wf.f, ", ", wf.interval, ")")

function Base.show(io::IO, mime::MIME"text/plain", wf::Waveform)
    xs = sample_clock(wf, 1e-3)
    plt = lineplot(
        xs, sample_values(wf);
        title=summary(wf),
        # TODO: decide the unit?
        xlabel="clock (μs)",
        ylabel="value (rad/µs)",
        compact=true,
    )
    return show(io, mime, plt)
end

function assert_interval_equal(lhs::Waveform, rhs::Waveform)
    lhs.interval == rhs.interval ||
        throw(ArgumentError("waveforms intervals are different cannot add them"))
end

function Base.:+(lhs::Waveform, rhs::Waveform)
    assert_interval_equal(lhs, rhs)
    return Waveform(lhs.interval) do t
        lhs.f(t) + rhs.f(t)
    end
end

function Base.:-(lhs::Waveform, rhs::Waveform)
    assert_interval_equal(lhs, rhs)
    return Waveform(lhs.interval) do t
        lhs.f(t) - rhs.f(t)
    end
end

function Base.:-(wf::Waveform)
    return Waveform(wf.interval) do t
        -wf.f(t)
    end
end

function append(wfs::Waveform...)
    last = wfs[1].interval.last
    checkpoints = Vector{typeof(last)}(undef, length(wfs))
    offsets = Vector{typeof(last)}(undef, length(wfs))

    checkpoints[1] = wfs[1].interval.last
    offsets[1] = 0
    @inbounds for idx in 2:length(wfs)
        last += span(wfs[idx].interval)
        checkpoints[idx] = last
        offsets[idx] = wfs[idx-1].interval.last - wfs[idx].interval.first
    end
    interval = Interval{Closed, Closed}(wfs[1].interval.first, last)

    return Waveform(interval) do t
        idx = 1
        while idx ≤ length(wfs) && t > checkpoints[idx]
            idx += 1
        end
        return wfs[idx].f(t - offsets[idx])
    end
end

function assert_clocks(clocks)
    issorted(clocks) || throw(ArgumentError("expect clocks to be sorted"))
    return
end

# this is for accessing the clocks and values
# in pulse smoothen, we may remove this if a more
# general version of the smoothen is implemented
struct PiecewiseLinear{T <: Real, Interp}
    clocks::Vector{T}
    values::Vector{T}
    interp::Interp

    function PiecewiseLinear(clocks::Vector{<:Real}, values::Vector{<:Real})
        assert_clocks(clocks)
        interp = LinearInterpolation(clocks, values)
        new{eltype(values), typeof(interp)}(clocks, values, interp)
    end
end

(f::PiecewiseLinear)(t::Real) = f.interp(t)

function piecewise_linear(;clocks::Vector{<:Real}, values::Vector{<:Real})
    return Waveform(PiecewiseLinear(clocks, values), first(clocks)..last(clocks))
end

function piecewise_constant(;
        clocks::Vector{<:Real}, values::Vector{<:Real},
        start::Real=first(clocks),
        stop::Real=last(clocks),
    )
    assert_clocks(clocks)
    return Waveform(start..stop) do t
        idx = findfirst(>(t), clocks) # we checked range
        isnothing(idx) && return values[end]
        return values[idx-1]
    end
end

function linear_ramp(;start::Real=0.0, stop::Real, start_value::Real, stop_value::Real)
    start ≤ stop || throw(ArgumentError("expect start ≤ stop, got start=$start, stop=$stop"))
    return Waveform(start..stop) do t
        (stop_value - start_value) / (stop - start) * (t - start) + start_value
    end
end

function constant(;start::Real=0.0, stop::Real, value::Real)
    return Waveform(start..stop) do t
        value
    end
end

function sinusoidal(;start::Real=0.0, stop::Real, amplitude::Real=zero(start))
    return Waveform(start..stop) do t
        amplitude * sin(t)
    end
end
