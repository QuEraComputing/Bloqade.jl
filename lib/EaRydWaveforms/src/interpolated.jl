"""
    InterpolatedWaveform{Interp <: Interpolations.AbstractExtrapolation} <: AbstractWaveform

Waveform defined by interpolations created by [Interpolations](https://github.com/JuliaMath/Interpolations.jl).
"""
struct InterpolatedWaveform{T <: Real, Interp <: Interpolations.AbstractExtrapolation{T}} <: AbstractWaveform
    interpolation::Interp
    duration::T
end

function InterpolatedWaveform(interpolation::Interpolations.AbstractExtrapolation)
    clocks = Interpolations.knots(interpolation)
    assert_interpolted_clocks(clocks)
    InterpolatedWaveform(interpolation, maximum(xs))
end

"""
    InterpolatedWaveform(clocks::Vector{T}, values::Vector{T}) where {T <: Real}

Create a linear interpolated waveform from `clocks` and `values`.
"""
function InterpolatedWaveform(clocks::Vector{T}, values::Vector{T}) where {T <: Real}
    assert_interpolted_clocks(clocks)
    interpolation = LinearInterpolation(clocks, values)
    return InterpolatedWaveform(interpolation, maximum(clocks))
end

Base.eltype(waveform::InterpolatedWaveform{T}) where T = T
duration(waveform::InterpolatedWaveform) = waveform.duration

function (waveform::InterpolatedWaveform)(t::Real, offset::Real=zero(t))
    assert_clock(t, duration(waveform), offset)
    return waveform.interpolation(t-offset)
end

function assert_interpolted_clocks(clocks::Vector)
    iszero(first(clocks)) ||
        throw(ArgumentError("clock must start from zero, got clocks[1] = $(clocks[1])"))
    return
end

# TODO: implement smoothening

# struct GaussianFilter
#     # parameters ?
# end

# function smoothen(waveform::InterpolatedWaveform, filter::FilterType=GausianFilter)
# end
