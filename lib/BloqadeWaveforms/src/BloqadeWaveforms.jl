module BloqadeWaveforms

using Intervals
using UnicodePlots
using LinearAlgebra
using Interpolations

export 
    Waveform,
    sample_values,
    sample_clock,
    piecewise_constant,
    piecewise_linear,
    linear_ramp,
    constant,
    sinusoidal,
    append,
    ..,
    # smooth
    smooth,
    Kernels

include("waveform.jl")
include("smooth.jl")

end
