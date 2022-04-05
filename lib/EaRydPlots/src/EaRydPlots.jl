module EaRydPlots

export draw, plt, bitstring_hist

using YaoArrayRegister
using YaoSubspaceArrayReg
using EaRydWaveforms
using PyPlot: plt

include("waveform.jl")
include("bitstring_hist.jl")
# using Makie
# using EaRydCore
# using EaRydWaveforms

# export bitstring_histgram, bitstring_histgram!, draw, draw!

# include("bitstring_hist.jl")
# include("waveform.jl")

end # EaRydPlots
