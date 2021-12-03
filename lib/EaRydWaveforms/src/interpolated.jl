# NOTE: the type in Interpolations is quite complicated to support all kinds of interpolation
# so instead of using FunctionWaveform directly let's just wrap a custom struct to simplify
# type parameters we are interested in (basically only `T` is useful)

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
    assert_interpolated_clocks(clocks)
    InterpolatedWaveform(interpolation, maximum(xs))
end

Base.eltype(waveform::InterpolatedWaveform{T}) where T = T
duration(waveform::InterpolatedWaveform) = waveform.duration

function (waveform::InterpolatedWaveform)(t::Real, offset::Real=zero(t))
    assert_clock(t, duration(waveform), offset)
    return waveform.interpolation(t-offset)
end

function assert_interpolated_clocks(clocks::Vector)
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
