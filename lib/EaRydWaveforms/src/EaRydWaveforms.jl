module EaRydWaveforms

using Intervals
using UnicodePlots
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
    ..

include("waveform.jl")

end
