using BloqadeLattices
using LinearAlgebra




abstract type AbstractBoundedLattice{D} <: AbstractLattice{D} end

struct GeneralLatticeTile{D,K,T,Tile <: AbstractTile} <: AbstractLattice{D}
    vectors::NTuple{D,NTuple{D,T}}
    sites::NTuple{K,NTuple{D,T}}
    tile::Tile
    PBC::Bool
end
function GeneralLatticeTile(vectors, sites, bounds)
    bounds = (bounds...,)
    vectors = (vectors...,)
    sites = (sites...,)

    @assert Dims(vectors)==Dims(bounds)
    

end


generate_sties(::GeneralLatticeTile) = ...

generate_rydberg_interactions(::GeneralLatticeTile) = ...
generate_rydberg_interactions(::AtomList) = ...

