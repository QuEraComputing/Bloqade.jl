module EaRydWaveforms

using Interpolations
using UnicodePlots

export InterpolatedWaveform,
    RampWaveform,
    SinusoidalWaveform,
    CompositeWaveform

include("waveforms.jl")
include("interpolated.jl")
include("ramp.jl")
include("sinusoidal.jl")
include("composite.jl")

end
