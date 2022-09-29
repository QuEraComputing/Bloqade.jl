abstract type AbstractSuperCell{D} end

struct Parallelepiped{D, T} <: AbstractSuperCell{D}
    # abstract matrix dimensions already guaranteed to be 2D
    bounds::AbstractMatrix{T}
    bounds_inv::AbstractMatrix{T}
end

function Parallelepiped(bounds::AbstractMatrix{T}) where {T}
    D = size(bounds,1)
    @assert size(bounds,1) == size(bounds,2) # ensure proper square matrix

    # return Parallelpiped
    return Parallelepiped{D,T}(bounds,inv(bounds))
end

# check if a point is in the tile
# Tile{D,T} is just concrete type based off the Supercell{D}
# NTuple gives has D elements of type T, representing point
in_range(x) = 0 ≤ x < 1 ? true : false
within_cell(cell::Parallelepiped{D, T},x::NTuple{D,T}) where {D,T} = all(in_range.(cell.bounds_inv * collect(x)))

function wrap_around(tile::Parallelepiped{D},x::NTuple{D,T}) where {D,T}

end

# distance(::Supercell{D},x::NTuple{D,T},y::NTuple{D,T}) = 

