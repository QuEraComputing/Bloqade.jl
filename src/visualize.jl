export vizconfig
using BloqadeLattices: AtomList

# TODO: generalize to color scale
function vizconfig(atoms::AtomList; config)
    colors = map(c -> iszero(c) ? "blue" : "red", config)
    img, (dx, dy) = BloqadeLattices.img_atoms(atoms; scale = 1.5, colors = colors, blockade_radius = 0)
    return img
end
