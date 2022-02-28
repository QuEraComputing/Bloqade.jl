module EaRydPlots

using Makie
using EaRydCore
using EaRydWaveforms

export bitstring_histgram, bitstring_histgram!, draw, draw!

include("bitstring_hist.jl")
include("waveform.jl")

end # EaRydPlots
