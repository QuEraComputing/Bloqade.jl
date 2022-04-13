using Intervals

"""
    Waveform{F, T <: Real}

Type for waveforms. `Waveform`s are defined
as a function combiend with a real number
duration.

# Fields

- `f`: a callable object.
- `duration`: a real number defines the duration of this waveform.
"""
struct Waveform{F, T <: Real}
    f::F
    duration::T
end

"""
    Waveform(f; duration::Real)

Create a `Waveform` object from callable `f`,
the unit of `duration` is `μs`.

# Example

```julia-repl
julia> Waveform(duration=1.5) do t
    2t+1
end
           ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀Waveform{_, Float64}⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀ 
           ┌────────────────────────────────────────┐ 
         4 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡠⠞⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
           │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⠞⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
           │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡠⠞⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
           │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⠞⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
           │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⠚⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
           │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⠞⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
           │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⠚⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
value (rad/µs) │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡠⠞⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
           │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⠞⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
           │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⠞⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
           │⠀⠀⠀⠀⠀⠀⠀⠀⣠⠞⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
           │⠀⠀⠀⠀⠀⠀⣠⠞⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
           │⠀⠀⠀⠀⣠⠞⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
           │⠀⠀⣠⠞⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
         1 │⣠⠞⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
           └────────────────────────────────────────┘ 
           ⠀0⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀clock (μs)⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀2⠀
```
"""
Waveform(f; duration::Real) = Waveform(f, duration)
Base.eltype(wf::Waveform) = typeof(wf.duration)

function Base.getindex(wf::Waveform, slice::Interval{<:Real, Closed, Closed})
    issubset(slice, 0..wf.duration) || throw(ArgumentError("slice is not in $(wf.duration) range, got $slice"))
    return Waveform(slice.last - slice.first) do t
        wf(t + slice.first)
    end
end

function (wf::Waveform)(t::Real, offset::Real=zero(t))
    t - offset < wf.duration || t - offset ≈ wf.duration || throw(ArgumentError(
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

function Base.show(io::IO, wf::Waveform)
    if get(io, :compact, false)
        print(io, "Waveform(", wf.f, ", ", wf.duration, ")")
    else
        print(io, "Waveform(_, ", wf.duration, ")")
    end
end

function Base.show(io::IO, mime::MIME"text/plain", wf::Waveform)
    clocks = sample_clock(wf)
    plt = lineplot(
        clocks, _rm_err.(sample_values(wf, clocks));
        title="Waveform{_, $(eltype(wf))}",
        # TODO: decide the unit?
        xlabel="clock (μs)",
        ylabel="value (rad/µs)",
        compact=true,
    )
    return show(io, mime, plt)
end

# NOTE: we don't plot error bar in terminal
function _rm_err(x)
    hasfield(typeof(x), :val) && return x.val
    return x
end

function assert_duration_equal(lhs::Waveform, rhs::Waveform)
    lhs.duration ≈ rhs.duration ||
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

function Base.:*(alpha::Number, wf::Waveform)
    return Waveform(wf.duration) do t
        alpha * wf.f(t)
    end
end

# let's assume they are communitive
function Base.:*(wf::Waveform, alpha::Number)
    return alpha * wf
end

Base.broadcastable(x::Waveform) = Ref(x)

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
    all(≥(0), clocks) || throw(ArgumentError("clocks must be non-nagative values"))
    iszero(first(clocks)) || throw(ArgumentError("the starting clock must be zero"))
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
        length(clocks) == length(values) || throw(ArgumentError("expect clocks has the same length as values"))
        interp = LinearInterpolation(clocks, values)
        new{eltype(values), typeof(interp)}(clocks, values, interp)
    end
end

(f::PiecewiseLinear)(t::Real) = f.interp(t)

struct PiecewiseConstant{T <: Real}
    clocks::Vector{T}
    values::Vector{T}

    function PiecewiseConstant(clocks::Vector{<:Real}, values::Vector{<:Real})
        assert_clocks(clocks)
        length(clocks) == length(values) + 1 || throw(ArgumentError("expect clocks has the same length as values"))
        new{eltype(values)}(clocks, values)
    end
end

function (f::PiecewiseConstant)(t::Real)
    idx = findfirst(>(t), f.clocks) # we checked range
    isnothing(idx) && return f.values[end]
    return f.values[idx-1]
end

"""
    piecewise_linear(;clocks, values)

Create a piecewise linear waveform.

# Keyword Arguments

- `clocks::Vector{<:Real}`: the list of clocks for the corresponding values.
- `values::Vector{<:Real}`: the list of values at each clock.

# Example

```julia-repl
julia> piecewise_linear(clocks=[0.0, 2.0, 3.0, 4.0], values=[0.0, 2.0, 2.0, 0.0])
                  ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀Waveform{_, Float64}⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀ 
                  ┌────────────────────────────────────────┐ 
                2 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡞⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⢧⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⠏⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⡄⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡴⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢹⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡞⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢇⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⠏⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⡆⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡴⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢱⡀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡞⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢇⠀⠀⠀⠀⠀│ 
   value (rad/µs) │⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⠏⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⡄⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⡴⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢹⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⢀⡞⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢧⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⣠⠏⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⡄⠀⠀│ 
                  │⠀⠀⠀⠀⡴⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢱⡀⠀│ 
                  │⠀⠀⢀⡞⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢇⠀│ 
                  │⠀⣠⠏⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⡆│ 
                0 │⡴⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢱│ 
                  └────────────────────────────────────────┘ 
                  ⠀0⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀clock (μs)⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀4⠀ 
```

"""
function piecewise_linear(;clocks::Vector{<:Real}, values::Vector{<:Real})
    iszero(first(clocks)) || throw(ArgumentError("the first clock time should be zero"))
    return Waveform(PiecewiseLinear(clocks, values), last(clocks))
end

"""
    piecewise_constant(;clocks, values, duration=last(clocks))

Create a piecewise constant waveform.

# Keyword Arguments

- `clocks::Vector{<:Real}`: the list of clocks for the corresponding values.
- `values::Vector{<:Real}`: the list of values at each clock.
- `duration::Real`: the duration of the entire waveform, default is the last clock.

# Example

```julia-repl
julia> piecewise_constant(clocks=[0.0, 0.2, 0.5], values=[0.0, 1.5, 3.1], duration=1.1)
                  ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀Waveform{_, Float64}⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀ 
                  ┌────────────────────────────────────────┐ 
                4 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⢰⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
   value (rad/µs) │⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⢠⠒⠒⠒⠒⠒⠚⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                0 │⣀⣀⣀⣸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  └────────────────────────────────────────┘ 
                  ⠀0⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀clock (μs)⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀2⠀ 
```
"""
function piecewise_constant(;
        clocks::Vector{<:Real}, values::Vector{<:Real},
        duration::Real=last(clocks),
    )
    assert_clocks(clocks)
    return Waveform(PiecewiseConstant(clocks, values), duration)
end

"""
    linear_ramp(;duration, start_value, stop_value)

Create a linear ramp waveform.

# Keyword Arguments

- `duration::Real`: duration of the whole waveform.
- `start_value::Real`: start value of the linear ramp.
- `stop_value::Real`: stop value of the linear ramp.

# Example

```julia-repl
julia> linear_ramp(;duration=2.2, start_value=0.0, stop_value=1.0)
                  ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀Waveform{_, Float64}⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀ 
                  ┌────────────────────────────────────────┐ 
                1 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⠞⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⠞⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⠞⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠞⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡴⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡴⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡴⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
   value (rad/µs) │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡴⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡴⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡴⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⠀⠀⢀⡴⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⠀⠀⢀⡴⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⠀⠀⢀⡴⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  │⠀⢀⡴⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                0 │⡴⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
                  └────────────────────────────────────────┘ 
                  ⠀0⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀clock (μs)⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀3⠀ 
```
"""
function linear_ramp(;duration::Real, start_value::Real, stop_value::Real)
    return Waveform(duration) do t
        (stop_value - start_value) / duration * t + start_value
    end
end

"""
    constant(;duration::Real, value::Real)

Create a constant waveform.

# Keyword Arguments

- `duration::Real`: duration of the whole waveform.
- `value::Real`: value of the constant waveform.
"""
function constant(;duration::Real, value::Real)
    return Waveform(duration) do t
        value
    end
end

"""
    sinusoidal(;duration::Real, amplitude::Real=one(start))

Create a sinusoidal waveform of the following expression.

```julia
amplitude * sin(t)
```

# Keyword Arguments

- `duration`: duration of the waveform.
- `amplitude`: amplitude of the sin waveform.
"""
function sinusoidal(;duration::Real, amplitude::Real=one(duration))
    return Waveform(duration) do t
        amplitude * sin(t)
    end
end
