module BloqadePlots

export draw, draw!, plt, bitstring_hist

using YaoArrayRegister
using YaoSubspaceArrayReg
using BloqadeWaveforms
using PythonCall

const plt = PythonCall.pynew()

function __init__()
    # copied from PyPlotCall.jl
    PythonCall.pycopy!(plt, pyimport("matplotlib.pyplot"))
end

include("waveform.jl")
include("bitstring_hist.jl")

end # BloqadePlots
