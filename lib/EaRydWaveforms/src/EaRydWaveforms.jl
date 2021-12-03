module EaRydWaveforms

using Intervals
using UnicodePlots
using Interpolations

export 
    ConstantWaveform,
    InterpolatedWaveform,
    PiecewiseConsant,
    PiecewiseLinear,
    RampWaveform,
    SinusoidalWaveform,
    CompositeWaveform,
    SlicedWaveform,
    SumWaveform,
    NegWaveform,
    duration,
    sample_values,
    ..

include("waveforms.jl")
include("interpolated.jl")
include("ramp.jl")
include("sinusoidal.jl")
include("composite.jl")
include("piecewise.jl")
include("sum.jl")
include("constant.jl")
include("slice.jl")

end
