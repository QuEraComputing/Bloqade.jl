struct FunctionWaveform{F, T <: Real} <: AbstractWaveform
    f::F
    duration::T
end

function FunctionWaveform(f; duration::Real)
    return FunctionWaveform(f, duration)
end

duration(wf::FunctionWaveform) = wf.duration
Base.eltype(::FunctionWaveform{<:Any, T}) where T = T

function (wf::FunctionWaveform)(t::Real, offset::Real=zero(t))
    assert_clock(t, duration(wf), offset)
    return wf.f(t - offset)
end

# this is for pretty printing & dispatch
struct Sinusoidal{T} <: Function
    amp::T
end

(f::Sinusoidal)(t::Real) = f.amp * sin(t)

struct Constant{T} <: Function
    amp::T
end

(f::Constant)(::Real) = f.amp

const SinusoidalWaveform{T} = FunctionWaveform{Sinusoidal{T}, T}
const ConstantWaveform{T} = FunctionWaveform{Constant{T}, T}

function SinusoidalWaveform(;duration::Real, amplitude::Real=1.0)
    return FunctionWaveform(Sinusoidal(amplitude), duration)
end

function ConstantWaveform(;duration::Real, amplitude::Real=1.0)
    return FunctionWaveform(Constant(amplitude), duration)
end
