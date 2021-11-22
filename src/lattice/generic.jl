struct GenericLattice
    rows
end

struct Square{T <: Real}
    width::T
    height::T
    spacing::T
end

function atom_positions(ltc::Square{Int}) # return list of atom positions
    list = []
    for x in 0:ltc.spacing:ltc.width, y in 0:ltc.spacing:ltc.height
        push!(list, RydAtom(x, y))
    end
    return list
end
