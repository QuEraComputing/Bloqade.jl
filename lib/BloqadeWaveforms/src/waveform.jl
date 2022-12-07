"""
    struct Waveform{F,T<:Real}

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


"""
    LinearAlgebra.norm(x::Waveform;p::Real=1)

Defines the norm function on [`Waveform`](@ref) type.
"""
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
            2π(2t+1)
        end
                    ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀Waveform{_, Float64}⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
                    ┌────────────────────────────────────────┐
                  4 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡠⠞⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                    │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⣠⠞⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                    │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡠⠞⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                    │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡠⠞⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
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

julia>
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
    
    # NOTE: t - offset is always smaller than wf.duration
    # but it may be slightly larger due to floating point error
    # so we need to always use the smaller of the two values
    return wf.f(min(t - offset, wf.duration))
end

"""
    sample_clock(wf::Waveform; offset::Real = zero(eltype(wf)), dt::Real = 1e-3)

Generates range of time values based on `wf`'s duration 
with `dt` time between each time value along with
`offset` time added to the beginning and end of the waveform's time span.

See also [`sample_values`](@ref)

```jldoctest; setup=:(using BloqadeWaveforms)
julia> wf = sinusoidal(duration=2, amplitude=2π*2.2); # create a waveform

julia> sample_clock(wf;) # range from 0.0 to 2.0 with step of 0.001 (default arg)
0.0:0.001:2.0

julia> sample_clock(wf; offset=0.1) # offset beginning and end by 0.1
0.1:0.001:2.1

julia> sample_clock(wf; dt = 2e-3) # set step size of 2e-3
0.0:0.002:2.0
```
"""
function sample_clock(wf::Waveform; offset::Real = zero(eltype(wf)), dt::Real = 1e-3)
    return offset:dt:wf.duration+offset
end

"""
    sample_values(wf::Waveform, clocks; offset::Real = zero(eltype(wf)))
    sample_values(wf::Waveform; offset::Real = zero(eltype(wf)), dt::Real = 1e-3)

Samples of waveform `wf` values obtainable by either providing an iterable `clocks`
containing exact time values to sample from or providing `offset` and `dt` values 
which specify the offset to add to the beginning and end of the waveforms time span and 
the step between time values. 

See also [`sample_clock`](@ref)

```jldoctest; setup:=using(using BloqadeWaveforms)
julia> wf = linear_ramp(duration=0.5, start_value=0.0, stop_value=2π*1.0);

julia> sample_values(wf,0.0:0.1:0.5) # sample waveform values from range
6-element Vector{Float64}:
 0.0
 1.2566370614359172
 2.5132741228718345
 3.7699111843077517
 5.026548245743669
 6.283185307179586

julia> sample_values(wf; dt=5e-2) #5e-2 time gap between each sampled valued
11-element Vector{Float64}:
 0.0
 0.6283185307179586
 1.2566370614359172
 1.8849555921538759
 2.5132741228718345
 3.141592653589793
 3.7699111843077517
 4.39822971502571
 5.026548245743669
 5.654866776461628
 6.283185307179586

```
"""
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
    issorted(clocks) || throw(ArgumentError("clocks must be sorted"))
    all(≥(0), clocks) || throw(ArgumentError("clocks must be non-negative values"))
    iszero(first(clocks)) || throw(ArgumentError("the starting clock value must be zero"))
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
        length(clocks) == length(values) || throw(ArgumentError("clocks must have the same length as values"))
        T = promote_type(eltype(clocks), eltype(values))
        clocks = Vector{T}(clocks); values = Vector{T}(values);
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
        length(clocks) == length(values) + 1 || throw(ArgumentError("clocks must have one more element than values"))
        T = promote_type(eltype(values), eltype(clocks))
        return new{T}(clocks, values)
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
julia> piecewise_linear(clocks=[0.0, 2.0, 3.0, 4.0], values=2π * [0.0, 2.0, 2.0, 0.0])
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
    iszero(first(clocks)) || throw(ArgumentError("the first clock element should be zero"))
    return Waveform(PiecewiseLinear(clocks, values), last(clocks))
end

"""
    piecewise_constant(;clocks, values, duration=last(clocks))

Create a piecewise constant waveform. 

The number of elements in `clocks` should be one greater than the number of elements in `values`. 
For example, if you wanted to define a piecewise constant waveform with the following:
* 0.0 2π ⋅ MHz from 0.0 to 0.5 μs
* 1.0 2π ⋅ MHz from 0.5 to 0.9 μs 
* 2.1 2π ⋅ MHz from 0.9 to 1.1 μs
It would be expressed as:
`piecewise_constant(clocks=[0.0, 0.5, 0.9, 2.1], values=[0.0, 1.0, 2.1])`

# Keyword Arguments

- `clocks::Vector{<:Real}`: the list of clocks for the corresponding values.
- `values::Vector{<:Real}`: the list of values at each clock.
- `duration::Real`: the duration of the entire waveform, default is the last clock.

# Example

```julia-repl
julia> piecewise_constant(clocks=[0.0, 0.2, 0.5, 0.9], values=2π * [0.0, 1.5, 3.1])
                    ┌────────────────────────────────────────┐
                  4 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                    │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                    │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                    │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡖⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒│
                    │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                    │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                    │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
   value (2π ⋅ MHz) │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                    │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                    │⠀⠀⠀⠀⠀⠀⠀⠀⢰⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠒⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                    │⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                    │⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                    │⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                    │⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                  0 │⣀⣀⣀⣀⣀⣀⣀⣀⣸⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                    └────────────────────────────────────────┘
                    ⠀0⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀clock (μs)⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀0.9⠀

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
julia> linear_ramp(;duration=2.2, start_value=0.0, stop_value=1.0 * 2π)
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

"""
    function (..)(first, last)

Exported from [`Intervals`](https://github.com/invenia/Intervals.jl), creates a closed interval from `first..last` 
and can be used with `Waveform` structs to obtain a slice of a Waveform's values, with the waveform slice's time adjusted to begin at 0 μs and
the duration being `last - first`.

# Example

```julia
julia> wf = Waveform(t->2.2*2π*sin(2π*t), duration = 2)
                    ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀Waveform{_, Int64}⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
                    ┌────────────────────────────────────────┐
                  3 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                    │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                    │⠀⠀⠀⡴⠋⠙⢦⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡴⠋⠙⢦⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                    │⠀⠀⡼⠁⠀⠀⠈⢧⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡼⠁⠀⠀⠈⢧⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                    │⠀⢰⠃⠀⠀⠀⠀⠈⡆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢰⠃⠀⠀⠀⠀⠈⡆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                    │⢀⡏⠀⠀⠀⠀⠀⠀⢸⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡏⠀⠀⠀⠀⠀⠀⢸⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                    │⡼⠀⠀⠀⠀⠀⠀⠀⠀⢧⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡼⠀⠀⠀⠀⠀⠀⠀⠀⢧⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
   value (2π ⋅ MHz) │⠧⠤⠤⠤⠤⠤⠤⠤⠤⠼⡦⠤⠤⠤⠤⠤⠤⠤⠤⢤⠧⠤⠤⠤⠤⠤⠤⠤⠤⠼⡦⠤⠤⠤⠤⠤⠤⠤⠤⢤│
                    │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢳⠀⠀⠀⠀⠀⠀⠀⠀⡞⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢳⠀⠀⠀⠀⠀⠀⠀⠀⡞│
                    │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⣇⠀⠀⠀⠀⠀⠀⢸⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⣇⠀⠀⠀⠀⠀⠀⢸⠁│
                    │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠸⡄⠀⠀⠀⠀⢀⠇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠸⡄⠀⠀⠀⠀⢀⠇⠀│
                    │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢳⡀⠀⠀⢀⡞⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢳⡀⠀⠀⢀⡞⠀⠀│
                    │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠳⣄⣠⠞⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠳⣄⣠⠞⠀⠀⠀│
                    │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠀⠀⠀⠀│
                 -3 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                    └────────────────────────────────────────┘
                    ⠀0⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀clock (μs)⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀2⠀

julia> wf[0.9..1.5]
                    ⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀Waveform{_, Float64}⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀
                    ┌────────────────────────────────────────┐
                  3 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                    │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                    │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣠⠤⠔⠒⠒⠒⠦⢤⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                    │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡤⠞⠉⠀⠀⠀⠀⠀⠀⠀⠀⠈⠙⠲⣄⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                    │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡴⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠑⢦⡀⠀⠀⠀⠀⠀⠀│
                    │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡴⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦⠀⠀⠀⠀⠀│
                    │⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⡴⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠳⣄⠀⠀⠀│
   value (2π ⋅ MHz) │⠀⠀⠀⠀⠀⠀⠀⠀⣠⠎⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠳⡀⠀│
                    │⠀⠀⠀⠀⠀⠀⢀⡜⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⢦│
                    │⠉⠉⠉⠉⠉⡽⠋⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉│
                    │⠀⠀⠀⣠⠞⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                    │⠀⢀⡴⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                    │⡴⠋⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                    │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                 -2 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
                    └────────────────────────────────────────┘
                    ⠀0⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀clock (μs)⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀0.6⠀

```

"""
function (..) end