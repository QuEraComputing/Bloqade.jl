# Shape

"a function of time: f(t::Real) describes the shape"
abstract type AbstractPulseShape end

struct Constant <: AbstractPulseShape
    value
end

struct Linear <: AbstractPulseShape
    start_value
    finish_value
end

struct JLFunctionCall <: AbstractPulseShape
    func
end

struct BinaryOp <: AbstractPulseShape
    head::Symbol
    lhs
    rhs
end

abstract type AbstractWaveform end

"""
an object add a duration to the shape
"""
struct Waveform <: AbstractWaveform
    shape
    duration
end

"""
compose duration together
"""
struct Piecewise <: AbstractWaveform
    waveforms::Vector{Waveform}
end

struct Smooth <: AbstractWaveform
    # must be piecewise(linear/constant)
    waveform::Vector{Piecewise}
    kernel::Symbol
    edge_pad_size::Int
end

# TODO: check this with Harry, if this should be named as "pulse"
# or it should be named as "channel", or something better

# attach the waveform 
struct Pulse # or WaveformList?
    waveforms::Vector{Any} # atom label => AbstractWaveform
    glob # global channel
end
