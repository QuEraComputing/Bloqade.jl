using Braket: AnalogHamiltonianSimulationTaskResult
using Colors, ColorSchemes

"""
    plot(object; kw...)

Plot a given object from Bloqade. The object
can be a lattice object, an [`AtomList`](@ref),
or [`Waveform`](@ref).

!!! tip

    This function is not exported to avoid conflicts
    with other plotting packages, it is recommended to
    always use it as `Bloqade.plot`.
"""
function plot end

function plot! end

export plot_densities

include("waveform.jl")
include("lattices.jl")
include("bitstring_hist.jl")
