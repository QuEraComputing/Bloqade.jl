module BloqadeWaveforms

using Unitful: Unitful, μs, MHz, Quantity
using Intervals
using UnicodePlots
using LinearAlgebra
using Interpolations
using QuadGK: quadgk
using Optim: optimize,Brent
using BloqadeExpr: default_unit
using OrderedCollections: OrderedDict

export Waveform,
    sample_values,
    sample_clock,
    piecewise_constant,
    piecewise_linear,
    linear_ramp,
    constant,
    sinusoidal,
    append,
    piecewise_linear_interpolate,
    piecewise_constant_interpolate,
    norm,
    ..,
    # smooth
    smooth,
    Kernels

include("waveform.jl")
include("smooth.jl")
include("interpolate.jl")

end
