module EaRydPlots

export draw, plt, bitstring_hist

using YaoArrayRegister
using YaoSubspaceArrayReg
using EaRydWaveforms
using PythonCall

const plt = PythonCall.pynew()

function __init__()
    # copied from PyPlotCall.jl
    PythonCall.pycopy!(plt, pyimport("matplotlib.pyplot"))
end

include("waveform.jl")
include("bitstring_hist.jl")
# using Makie
# using EaRydCore
# using EaRydWaveforms

# export bitstring_histgram, bitstring_histgram!, draw, draw!

# include("bitstring_hist.jl")
# include("waveform.jl")

end # EaRydPlots
