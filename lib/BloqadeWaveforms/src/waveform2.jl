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
    waveform::Piecewise
    kernel::Symbol
    edge_pad_size::Int
end

# TODO: check this with Harry, if this should be named as "pulse"
# or it should be named as "channel", or something better

struct Channel
    label # => atom
    # individual
    rabi_waveform
    detuning_waveform
    # address::UInt64 # rabi/detuning 0x10000123
    # waveforms # waveforms in the same channel will be execute in parallel
end

struct RabiPattern
    mask::Vector{UInt8}
    waveform
end

struct DetuningPattern
    amp::Vector{Float64}
    waveform
end

struct RydbergProgram
    rabi::Vector{RabiPattern}
    detuning::Vector{DetuningPattern}
end

prog = RydbergProgram()

quote
    Δ(t) = 2t
    Ω(t) = 3t

    rydberg_h(Δ, Ω)
end

Vector{DetuningPattern}

0x00010000
0x01110110
0x00010010

Δ(i,t) = Δ_i*δ(i)

wf = δ/10 *(0+0+....1) = δ/10
