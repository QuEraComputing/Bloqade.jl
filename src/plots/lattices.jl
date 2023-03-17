function plot(atoms::AtomList; kw...)
    return BloqadeLattices.img_atoms(atoms; kw...)
end

function plot_densities(atoms::AtomList, reg; color_scheme=ColorSchemes.bwr)
    plot(atoms, colors = "#" .* hex.(get(color_scheme, rydberg_density(reg), :clamp)))
end

function plot_densities(atoms::AtomList, task_res::AnalogHamiltonianSimulationTaskResult; color_scheme=ColorSchemes.bwr)
    plot(atoms, colors = "#" .* hex.(get(color_scheme, rydberg_density(task_res), :clamp)))
end