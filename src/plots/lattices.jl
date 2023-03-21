function plot(atoms::AtomList; kw...)
    return BloqadeLattices.img_atoms(atoms; kw...)
end

"""
    plot_densities(atoms::AtomList, densities::AbstractVector; color_scheme=ColorSchemes.bwr)

Given `atoms` which contains the coordinates for the atoms and a vector containing the Rydberg densities,
plot the Rydberg densities over the atoms.
"""
function plot_densities(atoms::AtomList, densities::AbstractVector; color_scheme=ColorSchemes.bwr)
    length(atoms) == length(densities) || throw(ArgumentError("The number of atoms is not equal to the number of Rydberg density values."))
    plot(atoms, colors = "#" .* hex.(get(color_scheme, densities, :clamp)))
end