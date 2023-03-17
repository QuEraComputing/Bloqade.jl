function plot(atoms::AtomList; kw...)
    return BloqadeLattices.img_atoms(atoms; kw...)
end

"""
    plot_densities(atoms::AtomList, reg; color_scheme = ColorSchemes.bwr)
    plot_densities(atoms::AtomList, task_res::AnalogHamiltonianSimulationTaskResult; color_scheme = ColorSchemes.bwr)

Given `atoms` which contains the coordinates for the atoms and either a register or `AnalogHamiltonianSimulationQuantumTaskResult` 
from Braket.jl, plot the Rydberg densities over the atoms.
"""
function plot_densities(atoms::AtomList, reg; color_scheme=ColorSchemes.bwr)
    plot(atoms, colors = "#" .* hex.(get(color_scheme, rydberg_density(reg), :clamp)))
end

function plot_densities(atoms::AtomList, task_res::AnalogHamiltonianSimulationTaskResult; color_scheme=ColorSchemes.bwr)
    plot(atoms, colors = "#" .* hex.(get(color_scheme, rydberg_density(task_res), :clamp)))
end