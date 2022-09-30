abstract type AbstractSuperCell{D} end

struct Parallelepiped{D, T} <: AbstractSuperCell{D}
    # abstract matrix dimensions already guaranteed to be 2D
    bounds::Matrix{T}
    bounds_inv::Matrix{T}
end

function Parallelepiped(bounds)
    D = size(bounds,1)
    bounds_inv = inv(bounds)
    # return Parallelpiped
    return Parallelepiped{D,eltype(bounds_inv)}(bounds,bounds_inv)
end

function Parallelepiped(bounds::T) where {T<:Real}
    bounds = fill(bounds,1,1)
    bounds_inv = inv(bounds)
    # return Parallelpiped
    return Parallelepiped{1,eltype(bounds_inv)}(bounds,bounds_inv)
end

# check if a point is in the tile
# Tile{D,T} is just concrete type based off the Supercell{D}
# NTuple gives has D elements of type T, representing point
in_range(x) = 0 ≤ x < 1 ? true : false
within_cell(cell::Parallelepiped{D, T},x::NTuple{D,T}) where {D,T} = all(in_range.(cell.bounds_inv * collect(x)))


distance(x,y) = √sum((x .- y) .^ 2)


function wrap_around(cell::Parallelepiped{D,T},x) where {D,T}
end


function distance(cell::Parallelepiped{D,T},x,y) where {D,T}
    x = wrap_around(cell,[x...,])
    y = wrap_around(cell,[y...,])
    dist = BloqadeLattices.distance(x,y)

    for a in 1:D
        shift = cell.bounds[:,a]
        for b in a:D
            shift .+= (b > a ? cell.bounds[:,b] : 0)
            shift_y = y .+ shift
            shift_x = x .+ shift
            dist = min(
                BloqadeLattices.distance(x,shift_y),
                BloqadeLattices.distance(shift_x,y),
                dist
            )
        end
    end

    return dist
end

