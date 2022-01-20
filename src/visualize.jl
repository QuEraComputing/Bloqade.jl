export vizconfig
using EaRydLattices: AtomList

# TODO: generalize to color scale
function vizconfig(atoms::AtomList; config)
    colors=map(c->iszero(c) ? "blue" : "red", config)
    img, (dx, dy) = EaRydLattices.img_atoms(atoms; scale=1.5, colors=colors, blockade_radius=0)
    img
end

