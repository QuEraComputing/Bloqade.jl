"""
    struct Waveform

Type for waveforms. `Waveform`s are defined
as a function combined with a real number
duration.

# Fields

- `f`: a callable object.
- `duration`: a real number defines the duration of this waveform; default unit is `μs`.
"""
struct Waveform{F,T<:Real}
    f::F
    duration::T

    function Waveform(f, duration)
        duration = default_unit(μs, duration)
        duration ≥ 0 || throw(ArgumentError("duration must be non-negative"))
        return new{typeof(f),typeof(duration)}(f, duration)
    end
end

function Base.:(==)(lhs::Waveform, rhs::Waveform)
    return lhs.duration == rhs.duration && lhs.f == rhs.f
end



function LinearAlgebra.norm(x::Waveform;p::Real=1)
    if isfinite(p) # p-norm 
        kernel = t->abs.(x(t)) .^ p
        area,area_error = quadgk(kernel,0,x.duration)
        return area^(1/p)
    elseif isinf(p) # infinite norm
        kernel = t -> -abs.(x(t))
        res = optimize(kernel,0,x.duration,Brent())
        return -minimum(res)
    else
        throw(OverflowError("argument ord must be finite real valued or infinity."))
    end
end

function Base.isapprox(lhs::Waveform,rhs::Waveform;atol::Real=0,rtol::Real = (atol>0 ? 0 : √eps()),p::Real=1)
    if lhs != rhs
            return LinearAlgebra.norm(lhs-rhs;p) < max(atol,rtol*max(LinearAlgebra.norm(lhs;p),LinearAlgebra.norm(rhs;p)))
    else
        return true
    end
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
value (2π ⋅ MHz) │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡠⠞⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
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

function Base.getindex(wf::Waveform, slice::Interval{<:Real,Closed,Closed})
    issubset(slice, 0 .. wf.duration) || throw(ArgumentError("slice is not in $(wf.duration) range, got $slice"))
    return Waveform(slice.last - slice.first) do t
        return wf(t + slice.first)
    end
end

function (wf::Waveform)(t::Real, offset::Real = zero(t))
    t - offset < wf.duration ||
        t - offset ≈ wf.duration ||
        throw(ArgumentError("t is not in range, expect $(offset) ≤ t ≤ $(wf.duration + offset), got $t"))
    return wf.f(t - offset)
end

function sample_clock(wf::Waveform; offset::Real = zero(eltype(wf)), dt::Real = 1e-3)
    return offset:dt:wf.duration+offset
end

function sample_values(wf::Waveform, clocks; offset::Real = zero(eltype(wf)))
    return [wf(t, offset) for t in clocks]
end

function sample_values(wf::Waveform; offset::Real = zero(eltype(wf)), dt::Real = 1e-3)
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
        clocks,
        _rm_err.(sample_values(wf, clocks) ./ (2π));
        title = "Waveform{_, $(eltype(wf))}",
        # TODO: decide the unit?
        xlabel = "clock (μs)",
        ylabel = "value (2π ⋅ MHz)",
        compact = true,
    )
    return show(io, mime, plt)
end

# NOTE: we don't plot error bar in terminal
function _rm_err(x)
    hasfield(typeof(x), :val) && return x.val
    return x
end

function assert_duration_equal(lhs::Waveform, rhs::Waveform)
    return lhs.duration ≈ rhs.duration || throw(ArgumentError("waveforms durations are different cannot add them"))
end

function Base.:+(lhs::Waveform, rhs::Waveform)
    assert_duration_equal(lhs, rhs)
    return Waveform(lhs.duration) do t
        return lhs.f(t) + rhs.f(t)
    end
end

function Base.:-(lhs::Waveform, rhs::Waveform)
    assert_duration_equal(lhs, rhs)
    return Waveform(lhs.duration) do t
        return lhs.f(t) - rhs.f(t)
    end
end

function Base.:-(wf::Waveform)
    return Waveform(wf.duration) do t
        return -wf.f(t)
    end
end

function Base.:*(alpha::Number, wf::Waveform)
    return Waveform(wf.duration) do t
        return alpha * wf.f(t)
    end
end

function Base.:/(alpha::Number, wf::Waveform)
    return Waveform(wf.duration) do t
        return alpha / wf.f(t)
    end
end

# let's assume they are communitive
function Base.:*(wf::Waveform, alpha::Number)
    return alpha * wf
end

function Base.:/(wf::Waveform, alpha::Number)
    return Waveform(wf.duration) do t
        return wf.f(t) / alpha
    end
end

Base.broadcastable(x::Waveform) = Ref(x)

"""
    append(wf::Waveform, wfs::Waveform...)

Append other waveforms to `wf` on time axis.
"""
function append(wf::Waveform, wfs::Waveform...)
    duration = wf.duration + sum(x -> x.duration, wfs)
    offsets = Vector{typeof(duration)}(undef, length(wfs))

    clock = wf.duration
    @inbounds for idx in eachindex(offsets)
        offsets[idx] = clock
        clock += wfs[idx].duration
    end

    return Waveform(duration) do t
        zero(wf.duration) ≤ t ≤ wf.duration && return wf(t)
        
        idx = searchsortedlast(offsets,t)
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
struct PiecewiseLinear{T<:Real,Interp}
    clocks::Vector{T}
    values::Vector{T}
    interp::Interp

    function PiecewiseLinear(clocks::Vector{<:Real}, values::Vector{<:Real})
        assert_clocks(clocks)
        length(clocks) == length(values) || throw(ArgumentError("expect clocks has the same length as values"))
        interp = LinearInterpolation(clocks, values)
        return new{eltype(values),typeof(interp)}(clocks, values, interp)
    end
end

function Base.:(==)(lhs::PiecewiseLinear, rhs::PiecewiseLinear)
    return lhs.clocks == rhs.clocks && lhs.values == rhs.values
end

function PiecewiseLinear(clocks::Vector{<:Quantity}, values::Vector{<:Quantity})
    return PiecewiseLinear(default_unit(μs, clocks), default_unit(MHz, values))
end

(f::PiecewiseLinear)(t::Real) = f.interp(t)

struct PiecewiseConstant{T<:Real}
    clocks::Vector{T}
    values::Vector{T}

    function PiecewiseConstant(clocks::Vector{<:Real}, values::Vector{<:Real})
        assert_clocks(clocks)
        length(clocks) == length(values) + 1 || throw(ArgumentError("expect clocks has one more element than values"))
        return new{eltype(values)}(clocks, values)
    end
end

function Base.:(==)(lhs::PiecewiseConstant, rhs::PiecewiseConstant)
    return lhs.clocks == rhs.clocks && lhs.values == rhs.values
end

function PiecewiseConstant(clocks::Vector{<:Quantity}, values::Vector{<:Quantity})
    return PiecewiseConstant(default_unit(μs, clocks), default_unit(MHz, values))
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
   value (2π ⋅ MHz) │⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⠏⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⡄⠀⠀⠀⠀│ 
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
function piecewise_linear(; clocks::Vector, values::Vector)
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
   value (2π ⋅ MHz) │⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
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
function piecewise_constant(; clocks::Vector, values::Vector, duration::Number = last(clocks))
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
   value (2π ⋅ MHz) │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡴⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
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
function linear_ramp(; duration, start_value, stop_value)
    duration = default_unit(μs, duration)
    start_value = default_unit(MHz, start_value)
    stop_value = default_unit(MHz, stop_value)

    return Waveform(duration) do t
        return (stop_value - start_value) / duration * t + start_value
    end
end

"""
    constant(;duration::Real, value::Real)

Create a constant waveform.

# Keyword Arguments

- `duration::Real`: duration of the whole waveform.
- `value::Real`: value of the constant waveform.
"""
function constant(; duration, value)
    duration = default_unit(μs, duration)
    value = default_unit(MHz, value)
    return Waveform(duration) do t
        return value
    end
end

"""
    sinusoidal(;duration::Real, amplitude::Real=one(start))

Create a sinusoidal waveform of the following expression.

```julia
amplitude * sin(2π*t)
```

# Keyword Arguments

- `duration`: duration of the waveform.
- `amplitude`: amplitude of the sin waveform.
"""
function sinusoidal(; duration, amplitude = one(duration))
    duration = default_unit(μs, duration)
    amplitude = default_unit(MHz, amplitude)
    return Waveform(duration) do t
        return amplitude * sin(2π * t)
    end
end
