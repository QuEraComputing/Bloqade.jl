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

struct Smooth <: AbstractPulseShape
    shape
    kernel::Symbol
    edge_pad_size::Int
end

"""
an object add a duration to the shape
"""
struct Duration
    shape
    duration
end

"""
compose duration together
"""
struct Piecewise
    waveforms::Vector{Duration}
end

# TODO: check this with Harry, if this should be named as "pulse"
# or it should be named as "channel", or something better
struct Pulse
    
end

struct Global
    waveforms
    glob_waveform
end

constant(1.0, duration=2.5)
linear([1.0, 2.0, 3.0], [2.0, 3.0, 1.0])
