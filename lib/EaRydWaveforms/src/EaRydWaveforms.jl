module EaRydWaveforms

using Interpolations

export PiecewiseLinearWaveform,
    RampWaveform,
    SinusoidalWaveform,
    CompositeWaveform

include("waveforms.jl")
include("interpolated.jl")
include("ramp.jl")
include("sinusoidal.jl")
include("composite.jl")

end
