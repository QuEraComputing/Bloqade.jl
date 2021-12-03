struct PiecewiseConsant{T <: Real} <: AbstractWaveform
    clocks::Vector{T}
    values::Vector{T}
    duration::T

    function PiecewiseConsant{T}(clocks, values, duration) where T
        assert_interpolated_clocks(clocks)
        new{T}(clocks, values, duration)
    end
end

# always use the precision of values
function PiecewiseConsant(clocks::Vector, values::Vector{T}, duration::Real = maximum(clocks)) where T
    return PiecewiseConsant{T}(clocks, values, duration)
end

# kwargs function
function PiecewiseConsant(;clocks::Vector, values::Vector, duration::Real = maximum(clocks)) where T
    return PiecewiseConsant(clocks, values, duration)
end

function (wf::PiecewiseConsant)(t::Real, offset::Real=zero(t))
    assert_clock(t, duration(wf), offset)
    idx = findfirst(>(t), wf.clocks) # we checked range
    isnothing(idx) && return wf.values[end]
    return wf.values[idx-1]
end

duration(wf::PiecewiseConsant) = wf.duration
Base.eltype(::PiecewiseConsant{T}) where T = T


struct PiecewiseLinear{T <: Real, Interp} <: AbstractWaveform
    interp_waveform::InterpolatedWaveform{T, Interp}

    PiecewiseLinear{T, Interp}(interp) where {T, Interp} = new{T, Interp}(interp) # raw constructor
    function PiecewiseLinear(clocks::Vector{Tc}, values::Vector{Tv}) where {Tc <: Real, Tv <: Real}
        assert_interpolated_clocks(clocks)
        interpolation = LinearInterpolation(clocks, values)
        return new{Tv, typeof(interpolation)}(InterpolatedWaveform(interpolation, maximum(clocks)))
    end
end

PiecewiseLinear(;clocks::Vector, values::Vector) = PiecewiseLinear(clocks, values)

duration(wf::PiecewiseLinear) = duration(wf.interp_waveform)
Base.eltype(wf::PiecewiseLinear) = eltype(wf.interp_waveform)
(wf::PiecewiseLinear)(t::Real, offset::Real=zero(t)) = wf.interp_waveform(t, offset)
