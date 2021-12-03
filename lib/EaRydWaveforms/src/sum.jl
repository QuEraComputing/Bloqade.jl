"""
    SumWaveform{Waveforms <: Tuple} <: AbstractWaveform
    SumWaveform(waveforms::AbstractWaveform...)

Summation of waveforms. This is the type that `+(waveform, waveform)`
returns, shouldn't be used directly as interface.
"""
struct SumWaveform{Waveforms <: Tuple} <: AbstractWaveform
    waveforms::Waveforms
    function SumWaveform(waveforms::Tuple)
        all(isequal(duration(first(waveforms))), map(duration, waveforms)) ||
            throw(ArgumentError("duration of the waveforms is different, cannot add them"))
        new{typeof(waveforms)}(waveforms)
    end
end

SumWaveform(waveforms::AbstractWaveform...) = SumWaveform(waveforms)

struct NegWaveform{Waveform <: AbstractWaveform} <: AbstractWaveform
    waveform::Waveform
end

duration(wf::SumWaveform) = duration(first(wf.waveforms))
duration(wf::NegWaveform) = duration(wf.waveform)

Base.eltype(wf::SumWaveform) = eltype(wf.waveforms[1])
Base.eltype(wf::NegWaveform) = eltype(wf.waveform)

function (wf::SumWaveform)(t::Real, offset::Real=zero(t))
    return sum(wf.waveforms) do waveform
        waveform(t, offset)
    end
end

function (wf::NegWaveform)(t::Real, offset::Real=zero(t))
    return -wf.waveform(t, offset)
end

Base.:+(lhs::AbstractWaveform, rhs::AbstractWaveform) = SumWaveform(lhs, rhs)

# always expand SumWaveform
Base.:+(lhs::SumWaveform, rhs::AbstractWaveform) = SumWaveform(lhs.waveforms..., rhs)
Base.:+(lhs::AbstractWaveform, rhs::SumWaveform) = SumWaveform(lhs, rhs.waveforms...)
Base.:+(lhs::SumWaveform, rhs::SumWaveform) = SumWaveform(lhs.waveforms..., rhs.waveforms...)

Base.:-(x::AbstractWaveform) = NegWaveform(x)
Base.:-(x::NegWaveform) = x.waveform
Base.:-(lhs::AbstractWaveform, rhs::AbstractWaveform) = lhs + (-rhs)
Base.:-(lhs::AbstractWaveform, rhs::NegWaveform) = lhs + rhs.waveform
Base.:-(lhs::NegWaveform, rhs::AbstractWaveform) = lhs + (-rhs) # disambiguity
Base.:-(lhs::NegWaveform, rhs::NegWaveform) = lhs + rhs.waveform
