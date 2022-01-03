using Intervals

struct Waveform{F, T <: Real}
    f::F
    duration::T
end

Waveform(f; duration::Real) = Waveform(f, duration)
Base.eltype(wf::Waveform) = typeof(wf.duration)

function Base.getindex(wf::Waveform, slice::Interval{<:Real, Closed, Closed})
    issubset(slice, 0..wf.duration) || throw(ArgumentError("slice is not in $(wf.duration) range, got $slice"))
    return Waveform(slice.last - slice.first) do t
        wf(t + slice.first)
    end
end

function (wf::Waveform)(t::Real, offset::Real=zero(t))
    t - offset ≤ wf.duration || throw(ArgumentError(
        "t is not in range, expect $(offset) ≤ t ≤ $(wf.duration + offset), got $t")
    )
    return wf.f(t - offset)
end

function sample_clock(wf::Waveform; offset::Real=zero(eltype(wf)), dt::Real=1e-3)
    return offset:dt:wf.duration+offset
end

function sample_values(wf::Waveform, clocks; offset::Real=zero(eltype(wf)))
    return [wf(t, offset) for t in clocks]
end

function sample_values(wf::Waveform; offset::Real=zero(eltype(wf)), dt::Real=1e-3)
    return sample_values(wf, sample_clock(wf; offset, dt))
end

Base.show(io::IO, wf::Waveform) = print(io, "Waveform(", wf.f, ", ", wf.duration, ")")

function Base.show(io::IO, mime::MIME"text/plain", wf::Waveform)
    clocks = sample_clock(wf)
    plt = lineplot(
        clocks, sample_values(wf, clocks);
        title="Waveform{_, $(eltype(wf))}",
        # TODO: decide the unit?
        xlabel="clock (μs)",
        ylabel="value (rad/µs)",
        compact=true,
    )
    return show(io, mime, plt)
end

function assert_duration_equal(lhs::Waveform, rhs::Waveform)
    lhs.duration == rhs.duration ||
        throw(ArgumentError("waveforms durations are different cannot add them"))
end

function Base.:+(lhs::Waveform, rhs::Waveform)
    assert_duration_equal(lhs, rhs)
    return Waveform(lhs.duration) do t
        lhs.f(t) + rhs.f(t)
    end
end

function Base.:-(lhs::Waveform, rhs::Waveform)
    assert_duration_equal(lhs, rhs)
    return Waveform(lhs.duration) do t
        lhs.f(t) - rhs.f(t)
    end
end

function Base.:-(wf::Waveform)
    return Waveform(wf.duration) do t
        -wf.f(t)
    end
end

"""
    append(wf::Waveform, wfs::Waveform...)

Append other waveforms to `wf` on time axis.
"""
function append(wf::Waveform, wfs::Waveform...)
    duration = wf.duration + sum(x->x.duration, wfs)
    offsets = Vector{typeof(duration)}(undef, length(wfs))

    clock = wf.duration
    @inbounds for (idx, wf) in enumerate(wfs)
        offsets[idx] = clock
        clock += wf.duration
    end

    return Waveform(duration) do t
        zero(wf.duration) ≤ t ≤ wf.duration && return wf(t)

        idx = 1
        while idx < length(wfs) && t > offsets[idx]
            idx += 1
        end
        return wfs[idx](t, offsets[idx])
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
    iszero(first(clocks)) || throw(ArgumentError("the first clock time should be zero"))
    return Waveform(PiecewiseLinear(clocks, values), last(clocks))
end

function piecewise_constant(;
        clocks::Vector{<:Real}, values::Vector{<:Real},
        duration::Real=last(clocks),
    )
    assert_clocks(clocks)
    return Waveform(duration) do t
        idx = findfirst(>(t), clocks) # we checked range
        isnothing(idx) && return values[end]
        return values[idx-1]
    end
end

function linear_ramp(;duration::Real, start_value::Real, stop_value::Real)
    return Waveform(duration) do t
        (stop_value - start_value) / duration * t + start_value
    end
end

function constant(;duration::Real, value::Real)
    return Waveform(duration) do t
        value
    end
end

function sinusoidal(;duration::Real, amplitude::Real=zero(start))
    return Waveform(duration) do t
        amplitude * sin(t)
    end
end
